# -*- coding: utf-8 -*-
"""
Created on 26 January 2025 14:57 (HSC DPhil students Room, Oxford)

@author: varad
"""
#all libraries that I will be using in the simulation

import numpy as np
import simpy 
import gradio as gr
import pandas as pd
import random
import csv
import matplotlib
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import simpy.resources
import os
import plotly.subplots as sp

#Non Modifiable variables
class parameters (object):
    '''
    This class contains all the constant non-modifiable parameters that will go into the model
    These mostly include service times and the time for which the simulation is supposed to run and halt etc
    '''

    experiment_no = 0 #incremented every time the main function is called
    number_of_runs = 3 #Total number of times the simulation will run for 1 experiment

    pt_per_day = 10 #Number of patients that visit the Gynae OPD every day (derived from AIIMS Bhopal Annual Report)(This could also be capped for a day if there are limited spots)
    pt_interarrival_time = 50 #Number of minutes in a working day / number of total patients expected during the day
    #modifiable factors, will be defined again in the relevant class
    #resources 
    #Staff
    #Since these are modifiable parameters, they are not implemented but only defined here, they are implemented to be inputted from the gradio app
    gynae_resident = None #Gynaecological residents that perform history and examination and pap smear collection in the routine OPDs in AIIMS Bhopal
    gynae_consulant = None #Gynaecological consultants that perform procedures such as LEEP (LEETZ) or hysterectomy etc
    pathologists = None #No of pathologists that interpret the pathology findings
    cytotechnicians = None #No of cytotechnicians that process the sample generated


    #Stuff
    num_pap_kits = 50 #number of kits that the hospital has to perform a pap smear (ayre's spatula, glass slide, preservative and box)
    num_pathology_consumables = 50 #consumables required for processing the pathological specimen
    num_colposcopy_consumables = 50 #consumables required for conducting colposcopy
    num_thermal_consumables = 50 #consumables required for LEEP
    num_ot_consumables = 50 #consumables required for hysterectomy


    #rooms (scheduled resource)
    num_colposcopy_rooms = 3 #number of colposcopy rooms
    num_ot_rooms = 3 #number of OT rooms

    #service times
    history_exam_time = 25 #time taken to complete history and examination per patient (Imp thing to remember here would be that this might change as the system adapts
                            #to an excess load)
    path_processing_time = 30 #time it takes from the sample is generated to the sample is prepared by cytotechnicians and ready for interpretation
    path_reporting_time = 30 # time taken by pathologists to report the results of a processed sample
    colposcopy_time = 25 #time taken to perform 1 colposcopy
    thermal_time = 25 #time taken for 1 Loop Electrosurgical Excision procedure
    hysterectomy_time = 50 #time taken for 1 hysterectomy

    #Epidemiological parameters
    screen_positivity_rate = 0.02 #% of positive samples (True positive + false positive / total samples)
    biopsy_rate = 0.4 # % of all colposcopies that undergo a biopsy
    biopsy_cin_rate = 0.6 # % of biopsies that are CIN
    biopsy_cacx_rate = 0.02 # % of biopsies that are CaCx
    follow_up_rate = 0.65 # % of women who follow up after a positive screen result (My own meta analysis + local data)

    run_time = 100000 #Time for which the entire simulation will run (in minutes)
class scheduled_resource(simpy.Resource):
    '''
    Extends the simpy.Resource object to include a resource that is only available during certain time of day and day of week
    '''
    def __init__(self, env, schedule, capacity):
        super().__init__(env, capacity)
        self.schedule = schedule # and integer list [0-6] for days of the week
        self.env = env

    def is_availeble (self):
        '''
        checks time of day and day of week and returns a boolean based on whether the resource is available at that time or not
        '''
        current_time = self.env.now
        week_minutes = 24 * 7 * 60 #minutes in a week
        day_minutes = 24 * 60 # minutes in a day

        current_day = int((current_time % week_minutes)/day_minutes) #first checks the number of minutes left in th  week then checks number of day of the week
        return current_day in self.schedule # returns a boolean whether the int current day is in schedule

    def request (self, *args, **kwargs):
        if self.is_availeble == False:
            self.env.process(self.wait_for_availability(*args, **kwargs))
        return super().request(*args, **kwargs)
    
    def wait_for_availability(self, *args, **kwargs):
        '''
        Creates a waiting process that waits for the resource to be available and then executes the request function
        '''
        while not self.is_availeble():
            #sees how much time is left for new day
            current_minutes = self.env.now
            day_minutes = 24 * 60
            minutes_till_next_day = day_minutes - (current_minutes/day_minutes)
            #wait for that much time
            yield self.env.timeout(minutes_till_next_day)
        #when it's the right time, execute the request
        request = super().request(*args, **kwargs)

        yield request
        return request


class ca_cx_patient (object):
    '''
    This class creates patients and declares their individual parameters that explains how they spent their time at the hospital
    These individual parameters will then be combined with others in the simulation to get overall estimates
    '''
    def __init__(self, pt_id):
        '''
        defines a patient and declares patient level variables to be recorded and written in a dataframe
        '''
        self.id = pt_id
        
        #declaring the variables to be recorded
        #putting them as zero to try an 
        self.time_at_entered = 0 #time when the patient entered into the OPD room
        self.time_at_screen_result = 0 #time when the patient first received the screening result
        self.time_at_colposcopy = 0 #time when the patient attended the colposcopy clinic
        self.time_at_treatment = 0 #time when patient got the treatment, either admission or surgery or LEEP or thermal/cryo
        self.time_at_exit = 0 #time when patient exits the system
        #need these values to calculate resource utilisation percentage
        self.history_examination_service_time = 0
        self.colposcopy_service_time = 0
        self.treatment_service_time = 0
        self.screen_sample_processing_time = 0
        self.screen_sample_reporting_time = 0
        self.biopsy_sample_processing_time = 0
        self.biopsy_sample_reporting_time = 0


        #need these values to calculate queue lengths 
        self.colposcopy_q_length =0
        self.treatment_q_length = 0
        self.screen_processing_q_length = 0
        self.screen_reporting_q_length = 0
        
        self.biopsy_processing_q_length = 0
        self.biopsy_reporting_q_length = 0


class Ca_Cx_pathway (object):
    '''
    This is the fake hospital. Defines all the processes that the patients will go through. Will record statistics for 1 simulation that will be later analyzed and clubbed with 
    results from 100 simulations
    '''
    def __init__(self, run_number, num_gynae_residents, num_gynae_consultants, num_pathologists, num_cytotechnicians, num_colposcopy_room = 3, num_ot_rooms = 2):
        self.env = simpy.Environment()

        #declaring number of modifiable resource capacity, non modifiable resources to be imported from the parameters class
        self.num_gynae_residents = num_gynae_residents
        self.num_gynae_consultants = num_gynae_consultants
        self.num_pathologists = num_pathologists
        self.num_cytotechnicians = num_cytotechnicians
        self.num_colposcopy_rooms = num_colposcopy_room
        self.num_ot_rooms = num_ot_rooms
        self.run_number = run_number

        self.colposcopy_schedule = [0,2,4] #list of integers form 0-6 for each day of the week that resource is available
        self.ot_schedule = [0,2,4]  #list of integers from 0-6 for each day of the week that resource is available
        
        self.pt_counter = 0 #acts as the UHID of the 0th patient
        self.run_number = self.run_number + 1
        #declaring resources
        #staff
        self.gynae_residents = simpy.Resource(self.env, capacity=num_gynae_residents)
        self.gynae_consultants = simpy.Resource(self.env, capacity=num_gynae_consultants)
        self.pathologist = simpy.Resource(self.env, capacity=num_pathologists)
        self.cytotechnician = simpy.Resource(self.env, capacity=num_cytotechnicians)

        #stuff
        self.pap_kit = simpy.Resource(self.env, capacity=parameters.num_pap_kits)
        self.pathology_consumables = simpy.Resource(self.env, capacity=parameters.num_pathology_consumables)
        self.colposcopy_consumables = simpy.Resource(self.env, capacity=parameters.num_colposcopy_consumables)
        self.thermal_consumables = simpy.Resource(self.env, capacity=parameters.num_thermal_consumables)
        self.ot_consumables = simpy.Resource(self.env, capacity=parameters.num_ot_consumables)

        #rooms (scheduled resource)
        self.colposcopy_room = scheduled_resource(self.env, self.colposcopy_schedule, capacity=parameters.num_colposcopy_rooms, )
        self.ot_room = scheduled_resource(self.env, self.ot_schedule, capacity = parameters.num_ot_rooms)
        
        #declaring a patient level dataframe to record patient KPIs - This is recorded at the individual level
        self.individual_results = pd.DataFrame({
            "UHID" : [],
            "Time_Entered_in System":[],
            "Screen_Processing_Q_Length" : [],
            "Screen_reporting_Q_Length" : [],
            "Time_at_screening_result":[],
            "Colposcopy_Q_Length" :[],
            "Time_at_colposcopy" : [],
            "Biopsy_Processing_Q_Length" : [],
            "Biopsy_Reporting_Q_Length" : [],
            "Treatment_Q_length" :[],
            "Time_at_treatment" : [],
            "History_and_Examination_time": [], #also recording service times as they will ultimately be added up to calculate resource utilisation percentage
            "Screen_processing_time":[],
            "Screen_reporting_time":[],
            "Biopsy_processing_time":[],
            "Biopsy_reporting_time":[],
            "Colposcopy_time":[],
            "Treatment_time":[],
            "Exit_time":[]
        })

        #Declaring individual results processing variables
        #time intervals between important points
        self.time_to_screen_result = 0 #during analysis, need to only consider those patients who actually did undergo these procedures
        self.time_to_colposcopy = 0    #as not all patients will undergo all the processes. might include some drop na function, but shouldn't be too much of a problem
        self.time_to_treatment = 0
        self.total_time_in_system = 0
        

        
        #declaring system KPIs to be measured at the run level.
        #Queue lengths for different processes
        self.max_q_len_screen_processing = 0
        self.max_q_len_screen_reporting = 0
        self.max_q_len_colposcopy = 0
        self.max_q_len_biopsy_processing = 0
        self.max_q_len_biopsy_reporting = 0
        self.max_q_len_treatment = 0

        #Resource utilization percentages
        self.gynae_residents_utilisation = 0 #by adding service times of all the processes where these resources are required. 
        self.gynae_consultants_utlisation = 0
        self.cytotechnician_utilisation = 0
        self.pathologist_utilisation = 0

        

    def is_within_working_hours(self):
        '''
        checks whether the current simulation time is within working hours and returns a boolean
        '''
        current_sim_mins = self.env.now
        day_mins = 24*60
        current_sim_hour = int((current_sim_mins%day_mins)/60)
        return 8 < current_sim_hour < 17


    
    def gen_patient_arrival(self):
        '''
        Generates a fictional patient according to a distribution, they undergo and OPD, this generates a sample which undergoes processing, after results are
        conveyed, if positive, patient only then moves on to the next step i.e. colposcopy.
        '''
        while True:
        #check time of day, 
            #if self.is_within_working_hours:
        #if time of day is appropriate then generate the patient
                
            self.pt_counter += 1
            screening_patient = ca_cx_patient(self.pt_counter)
                
                
                #print("Patient generates", self.patient.pt_id)
                #here we will need to generate all the samples for the patient, even if they don't get created later on
                #reason for that is it's okay if the reading is 0 or NaN, the code just won't work if there is no object to begin with

        #record necessary timepoints
            screening_patient.time_at_entered = self.env.now
        #patient moves to the OPD
            self.env.process(self.history_examination(screening_patient))
        #time for next patient arrival
            wait_time_for_next_pt = random.expovariate(1/parameters.pt_interarrival_time)
            yield self.env.timeout(wait_time_for_next_pt)

    def history_examination(self, patient):
        '''
        Patient undergoes history and examination and in the process also generates the screening sample
        '''
        #request for a resident and consumables for sample collection and wait for them to be available
        with self.gynae_residents.request() as gynae_res, self.pap_kit.request() as pap, self.pathology_consumables.request() as path_consum :
            yield gynae_res and pap and path_consum
        
        #patient undergoes history, examination and sample collection
            history_examination_time = random.triangular(parameters.history_exam_time/2, parameters.history_exam_time, parameters.history_exam_time *2 )
            patient.history_examination_service_time = history_examination_time
            yield self.env.timeout(history_examination_time)

            #New implementation different than the previous one
            #The sample goes to processing and reporting function which generates a boolean which decides whether the patient moves on or not
            screening_sample_gen = self.env.process(self.screening(patient))
            
            screen_result = yield screening_sample_gen
            
            
            if screen_result:
                self.env.process(self.call_for_follow_up(patient))
            else:
                #patient exits the system
                patient.time_at_exit = self.env.now
                self.add_to_individual_results(( patient))
            
            #generate a screening sample
            #self.pt_screening_sample = screen_sample(self.patient.pt_id) #screen sample id is the same as the patient id
            #print("Screen Sample generated", self.pt_screening_sample.screen_sample_id)
            #here we will need to generate all the samples for the patient, even if they don't get created later on
            #reason for that is it's okay if the reading is 0 or NaN, the code just won't work if there is no object to begin with
            #self.pt_biopsy_sample = biopsy_sample(self.patient.pt_id) #generate a biopsy sample that will go for processing  
            #print("biopsy sample generated", self.pt_biopsy_sample.biopsy_sample_id)
        
            #sample goes on for processing
            #self.env.process(self.screen_sample_processing())

    def screening (self, patient):
        '''
        This function simulation the processing and reporting of screen samples and returns a boolean whether the result is positive or negative
        '''
        patient.screen_processing_q_length = len(self.cytotechnician.queue)
        
        with self.cytotechnician.request() as cytotec, self.pathology_consumables.request() as scr_proc_consum:
            yield cytotec and scr_proc_consum

            screen_sample_processing_time = random.triangular(parameters.path_processing_time/2, parameters.path_processing_time, parameters.path_processing_time *2)
            patient.screen_sample_processing_time = screen_sample_processing_time
            yield self.env.timeout(screen_sample_processing_time)
            
            
            patient.screen_reporting_q_length = len(self.pathologist.queue)
            with self.pathologist.request() as path:
                yield path
                screen_sample_reporting_time = random.triangular(parameters.path_reporting_time/2, parameters.path_reporting_time, parameters.path_reporting_time * 2)
                patient.screen_sample_reporting_time = screen_sample_reporting_time #record this for resource utilisation %
                yield self.env.timeout(screen_sample_reporting_time)
                
                patient.time_at_screen_result = self.env.now

                if random.random() < parameters.screen_positivity_rate:
                    return True #if sample is positive
                else: 
                    return False # if sample is negative

    def screen_sample_processing(self):
        '''
        Sample undergoes processing
        '''
        #queue length for processing
        self.pt_screening_sample.screen_processing_q_length = len(self.cytotechnician.queue)
        #request resources and wait for them to be available
        with self.cytotechnician.request() as cytotec, self.pathology_consumables.request() as scr_proc_consum:
            yield cytotec and scr_proc_consum
        
        #sample undergoes processing
            screen_sample_processing_time = random.triangular(parameters.path_processing_time/2, parameters.path_processing_time, parameters.path_processing_time *2)
            self.pt_screening_sample.screen_sample_processing_time = screen_sample_processing_time
            yield self.env.timeout(screen_sample_processing_time)
        #sample goes for reporting
            self.env.process(self.screen_sample_reporting())

    def screen_sample_reporting(self):
        '''
        Processed sample is interpreted and reported by pathologist
        '''
        #measure queue length for every patient that comes (max of this column will be the max queue length)
        self.pt_screening_sample.screen_reporting_q_length = len(self.pathologist.queue)
        #request for a pathologist and wait until 
        with self.pathologist.request() as path:
            yield path

        #record the current time as an important milestone
            self.patient.time_at_screen_result = self.env.now
        
        #sample undergoes reporting
            screen_sample_reporting_time = random.triangular(parameters.path_reporting_time/2, parameters.path_reporting_time, parameters.path_reporting_time * 2)
            self.pt_screening_sample.screen_sample_reporting_time = screen_sample_reporting_time #record this for resource utilisation %
            yield self.env.timeout(screen_sample_reporting_time)

        #if sample is positive, move on to follow up, otherwise terminate
            if random.random() < parameters.screen_positivity_rate:
                self.env.process(self.call_for_follow_up())
            else:
                #patient exits the system
                self.patient.time_at_exit = self.env.now
                #add data to the df
                Ca_Cx_pathway.add_to_individual_results(self)

    def call_for_follow_up (self, patient):
        '''
        Gynaecology residents 
        '''
        #no waiting time for this as it is quite instant.
        #request a gynae_res (later on could modify to include a receptionist or another health cadre)
        with self.gynae_residents.request() as gynae_res:
            yield gynae_res
        
        # whether the patient returns or not
            if random.random() < parameters.follow_up_rate:
            #patient goes on for colposcopy
                self.env.process(self.colposcopy(patient))
        #instantaneous process so no timeout really and also not a service
            else:
                #patient exits the system
                patient.time_at_exit = self.env.now
                #add to df
                self.add_to_individual_results(patient)
            
    def colposcopy(self, patient):
        '''
        Patient that was generated undergoes colposcopy
        '''
        #here, the entity requests two different resources, it's waiting time or queue length will be decided by whatever is less available. 
        # 1 small caveat here is that service times for different resources are different, so a larger queue doesn't necessarily mean a longer waiting time
        # we're not measuring waiting time but only time between events as that is a much more relevant indicator for implementation decisions.
        colpo_q_len_list = [len(self.gynae_consultants.queue), len(self.colposcopy_room.queue) ]
        patient.colposcopy_q_length = max(colpo_q_len_list)
        
        #requests for a consultant, consumables and a room
        with self.gynae_consultants.request() as gynae_consul, self.colposcopy_consumables.request() as gynae_consumables, self.colposcopy_room.request() as colpo_room:
            yield gynae_consul and gynae_consumables and colpo_room

        #Record time at colposcopy
            patient.time_at_colposcopy = self.env.now

        #patient undergoes colposcopy
            colposcopy_service_time = random.triangular(parameters.colposcopy_time/2, parameters.colposcopy_time, parameters.colposcopy_time *2)
            patient.colposcopy_service_time = colposcopy_service_time
            yield self.env.timeout(colposcopy_service_time)

            

            biopsy_sample_gen = self.env.process(self.biopsy(patient))

            biopsy_result = yield biopsy_sample_gen

            if biopsy_result == 1:
                self.env.process(self.thermal_ablation(patient))
            elif biopsy_result == 2:
                self.env.process(self.hysterectomy(patient))
            else:
                patient.time_at_exit = self.env.now
                self.add_to_individual_results((patient))

    def biopsy(self, patient):
        '''
        implementation is very similar to the screening function
        '''
        patient.biopsy_processing_q_length = len(self.cytotechnician.queue)
        with self.cytotechnician.request() as cytotec, self.pathology_consumables.request() as path_consum:
            yield cytotec and path_consum

            biopsy_sample_processing_time = random.triangular(parameters.path_processing_time/2, parameters.path_processing_time, parameters.path_processing_time * 2)    
            patient.biopsy_sample_processing_time = biopsy_sample_processing_time
            yield self.env.timeout(biopsy_sample_processing_time)

            patient.biopsy_reporting_q_length = len(self.pathologist.queue)
            with self.pathologist.request() as path:
                yield path

                biopsy_sample_reporting_time = random.triangular(parameters.path_reporting_time/2, parameters.path_reporting_time, parameters.path_reporting_time *2)
                patient.biopsy_sample_reporting_time = biopsy_sample_reporting_time
                yield self.env.timeout(biopsy_sample_reporting_time)

                if random.random() < parameters.biopsy_cin_rate:
                    return 1
                elif parameters.biopsy_cin_rate < random.random() < parameters.biopsy_cacx_rate:
                    return 2
                else:
                    return 3
                    
    def biopsy_sample_processing(self):
        '''
        Biopsy sample if prepared undergoes processing
        '''
        #queue length for processing
        self.pt_biopsy_sample.biopsy_processing_q_length = len(self.cytotechnician.queue)
        #requests a cytotechnicians and consumables
        with self.cytotechnician.request() as cytotec, self.pathology_consumables.request() as path_consum:
            yield cytotec and path_consum

            #biopsy sample undergoes processing
            biopsy_sample_processing_time = random.triangular(parameters.path_processing_time/2, parameters.path_processing_time, parameters.path_processing_time * 2)    
            self.pt_biopsy_sample.biopsy_sample_processing_time = biopsy_sample_processing_time
            yield self.env.timeout(biopsy_sample_processing_time)

            #biopsy sample goes for reporting
            self.env.process(self.biopsy_sample_reporting())

    def biopsy_sample_reporting(self):
        '''
        Biopsy sample if taken undergoes reporting after processing
        '''
        #queue length for reporting
        self.pt_biopsy_sample.biopsy_reporting_q_length = len(self.pathologist.queue)
        #requests a pathologist
        with self.pathologist.request() as path:
            yield path

            #biopsy sample undergoes reporting
            biopsy_sample_reporting_time = random.triangular(parameters.path_reporting_time/2, parameters.path_reporting_time, parameters.path_reporting_time *2)
            self.pt_biopsy_sample.biopsy_sample_reporting_time = biopsy_sample_reporting_time
            yield self.env.timeout(biopsy_sample_reporting_time)

            #depending on the diagnosis, patient either goes for thermal ablation or hysterectomy (currently, only making 2 options available, have the option of adding more on later)
            biopsy_result = random.random()
            if biopsy_result < parameters.biopsy_cin_rate:
                self.env.process(self.thermal_ablation()) #diagnosed with CIN
            
            elif parameters.biopsy_cin_rate < biopsy_result < parameters.biopsy_cacx_rate:
                self.env.process(self.hysterectomy()) #diagnosed with cervical cancer
            
            else:
                self.patient.time_at_exit = self.env.now #patient exits the system
                #add data to the df
                Ca_Cx_pathway.add_to_individual_results(self)

    def thermal_ablation(self, patient):
        '''
        If indicated, pt undergoes thermal ablation
        '''
        thermal_q_len_list = [len(self.gynae_consultants.queue), len(self.colposcopy_room.queue) ]
        patient.treatment_q_length = max(thermal_q_len_list)
        
        #requests resources required for thermal ablation
        with self.gynae_consultants.request() as gynae_consul, self.thermal_consumables.request() as thermal_consum, self.colposcopy_room.request() as colpo_room:
            yield gynae_consul and thermal_consum and colpo_room

            patient.time_at_treatment = self.env.now
        #patient undergoes thermal ablation
            thermal_ablation_time = random.triangular(parameters.thermal_time/2, parameters.thermal_time, parameters.thermal_time *2)
            patient.treatment_service_time = thermal_ablation_time
            yield self.env.timeout(thermal_ablation_time)

        #patient exits the system
            
            patient.time_at_exit = self.env.now
            #add to df
            self.add_to_individual_results(patient)

    def leep (self):
        '''
        if indicated, patient undergoes LEEP
        '''
        #Not being implemented in this first version of the model
        pass 
    def hysterectomy (self, patient):
        '''
        if indicated, patient undergoes hysterectomy
        '''
        hyst_q_len_list = [len(self.gynae_consultants.queue), len(self.ot_room.queue)]
        patient.treatment_q_length = max(hyst_q_len_list)
        
        #request for a ot room and other equipment
        with self.gynae_consultants.request() as gynae_consul, self.ot_consumables.request() as ot_consum, self.ot_room as ot_room:
            yield gynae_consul and ot_consum and ot_room

            self.patient.time_at_treatment = self.env.now    
            #patient undergoes surgery
            hysterectomy_time = random.triangular(parameters.hysterectomy_time/2, parameters.hysterectomy_time, parameters.hysterectomy_time *2)
            self.patient.treatment_service_time = hysterectomy_time
            yield self.env.timeout(hysterectomy_time)

            #patient exits the system
            
            patient.time_at_exit = self.env.now
            #adding everything to the dataframe
            self.add_to_individual_results(patient)

    
    def add_to_individual_results (self, patient):
        '''
        To add a row to a df, we need to pass an argument that adds in all 10-12 columns together even if we want to add just one cell
        Hence to make my job easier, writing a function that does this in every function without having to write too much.
        '''
        df_to_add = pd.DataFrame({
            "UHID" : [patient.id],
            "Time_Entered_in System":[patient.time_at_entered],
            "Screen_Processing_Q_Length" : [patient.screen_processing_q_length],
            "Screen_reporting_Q_Length" : [patient.screen_reporting_q_length],
            "Time_at_screening_result":[patient.time_at_screen_result],
            "Colposcopy_Q_Length":[patient.colposcopy_q_length],
            "Time_at_colposcopy" : [patient.time_at_colposcopy],
            "Biopsy_Processing_Q_Length" : [patient.biopsy_processing_q_length],
            "Biopsy_Reporting_Q_Length" : [patient.biopsy_reporting_q_length],
            "Treatment_Q_length":[patient.treatment_q_length],
            "Time_at_treatment" : [patient.time_at_treatment],
            "History_and_Examination_time": [patient.history_examination_service_time], #also recording service times as they will ultimately be added up to calculate resource utilisation percentage
            "Screen_processing_time":[patient.screen_sample_processing_time],
            "Screen_reporting_time":[patient.screen_sample_reporting_time],
            "Biopsy_processing_time":[patient.biopsy_sample_processing_time],
            "Biopsy_reporting_time":[patient.biopsy_sample_reporting_time],
            "Colposcopy_time":[patient.colposcopy_service_time],
            "Treatment_time":[patient.treatment_service_time],
            "Exit_time":[patient.time_at_exit]
            
        })
        df_to_add.set_index('UHID', inplace= True)
        self.individual_results = pd.concat([self.individual_results, df_to_add]) #throws syntax error that I should not use the _ sign, we'll see
    

    def individual_results_processor(self):
        '''
        Processes the individual results dataframe by adding columns from which KPI's can be calculated
        '''
        #Calculating time between important events
        self.individual_results['Time_to_screen_results'] = self.individual_results['Time_at_screening_result'] - self.individual_results['Time_Entered_in System']
        self.individual_results['Time_to_Colposcopy'] = self.individual_results['Time_at_colposcopy'] - self.individual_results['Time_Entered_in System']
        self.individual_results['Time_to_Treatment'] = self.individual_results['Time_at_treatment'] - self.individual_results['Time_Entered_in System']
        self.individual_results['Total_time_in_system'] = self.individual_results['Exit_time'] - self.individual_results['Time_Entered_in System']
        
        #Calculating service times for different resources
        
        self.individual_results['Gynae_res_busy_time'] = self.individual_results['History_and_Examination_time']
        self.individual_results['Cytotech_busy_time'] = self.individual_results['Screen_processing_time'] + self.individual_results['Biopsy_processing_time']
        self.individual_results['Pathologist_busy_time'] = self.individual_results['Screen_reporting_time'] + self.individual_results['Biopsy_reporting_time']
        self.individual_results['Gynae_consul_busy_time'] = self.individual_results['Colposcopy_time'] + self.individual_results['Treatment_time']

    
            
    def KPI_calculator(self):
        '''
        Function that calculates the various KPIs from an individual run from the different columns of the individual results dataframe
        These are KPIs for a signle run
        '''
        #max q lengths
        self.max_q_len_screen_processing = self.individual_results['Screen_Processing_Q_Length'].max()
        self.max_q_len_screen_reporting = self.individual_results['Screen_reporting_Q_Length'].max()
        self.max_q_len_colposcopy = self.individual_results['Colposcopy_Q_Length'].max()
        self.max_q_len_biopsy_processing = self.individual_results['Biopsy_Processing_Q_Length'].max()
        self.max_q_len_biopsy_reporting = self.individual_results['Biopsy_Reporting_Q_Length'].max()
        self.max_q_len_treatment = self.individual_results['Treatment_Q_length'].max()
        
        #resource utilisation percentages
        self.gynae_residents_utilisation = self.individual_results['Gynae_res_busy_time'].sum()/(parameters.run_time * self.num_gynae_residents)
        self.cytotechnician_utilisation = self.individual_results['Cytotech_busy_time'].sum()/(parameters.run_time * self.num_cytotechnicians)
        self.gynae_consultants_utlisation = self.individual_results['Gynae_consul_busy_time'].sum() / (parameters.run_time * self.num_gynae_consultants)
        self.pathologist_utilisation = self.individual_results['Pathologist_busy_time'].sum() / (parameters.run_time * self.num_pathologists)

        #median time to important events
        #creating temp df and dropping rows with negative values for specific columns
        temp_colpo_time_df = self.individual_results['Time_to_Colposcopy'][self.individual_results['Time_to_Colposcopy'] >0]
        temp_treatmet_time_df = self.individual_results['Time_to_Treatment'][self.individual_results['Time_to_Treatment'] >0]
        #now putting the median method onto that limited dataset
        self.med_time_to_scr_res = self.individual_results['Time_to_screen_results'].median()
        self.med_time_to_colpo = temp_colpo_time_df.median()
        self.med_time_to_treatment = temp_treatmet_time_df.median()
        self.med_tot_time_in_system = self.individual_results['Total_time_in_system'].median()


    def export_row_to_csv(self):
        '''
        Creates a new dataframe with trial results and exports a single row to that dataframe after each run
        '''
        with open ('kpi_trial_results.csv', 'a')as f:
            writer = csv.writer(f, delimiter= ',')
            row_to_add = [
                self.run_number,
                self.max_q_len_screen_processing,
                self.max_q_len_screen_reporting,
                self.max_q_len_colposcopy,
                self.max_q_len_biopsy_processing,
                self.max_q_len_biopsy_reporting,
                self.max_q_len_treatment,

                self.gynae_residents_utilisation,
                self.gynae_consultants_utlisation,
                self.pathologist_utilisation,
                self.cytotechnician_utilisation,

                self.med_time_to_scr_res,
                self.med_time_to_colpo,
                self.med_time_to_treatment,
                self.med_tot_time_in_system
            ]
            writer.writerow(row_to_add)

    def run(self):
        '''
        Runs the simulation and calls the generator function.
        '''
        self.env.process(self.gen_patient_arrival())
        self.env.run(until= parameters.run_time)
        #print(self.individual_results)
        self.individual_results_processor()
        self.individual_results.to_csv('individual_results.csv')
        self.KPI_calculator()
        self.export_row_to_csv()


class summary_statistics(object):
    '''
    This class will define methods that will calculate aggregate statistics from 100 simulations and append the results onto a new spreadsheet which will be used to append results
    from 100 simulations for different number of independent variables (such as patients)
    '''
    def __init__(self):
        pass
    
    
          
    def gen_final_summary_table (self):
        '''
        Generates a table, essentially a row of summary statistics for 100 runs with a particular initial setting.
        '''
        with open ('final_summary_table.csv', 'w') as f:
            writer = csv.writer(f, delimiter= ',')
            column_headers = [
                "Experiment_No" ,
                "Max_Scr_Proc_Q_len" ,
                'Max_Scr_Rep_Q_len' ,
                'Max_Colpo_Q_len' ,
                "Max_Biop_Proc_Q_Len" ,
                "Max_Biop_Rep_Q_Len" ,
                "Max_T/t_Q_len" ,

                #resource utilisation %
                'Gynae_Res_%_util',
                'Gynae_consul_%_util',
                'Path_%_util',
                'Cytotec_%_util',

                #Time between important events
                'Time_to_screening_results',
                'Time_to_colposcopy',
                'Time_to_treatment',
                'Total_time_in_system' ]
            
            writer.writerow(column_headers)
        
    def calculate_summary_statistics(self):
        '''
        Calculates summary statistic from 100 runs (or whatever the number of runs is specified) from the kpi_trial_results table which will then later on be added
        onto the final_summary_table csv
        '''
        filepath = 'kpi_trial_results.csv'
        df_to_read = pd.read_csv(filepath)
        self.max_scr_proc_q_len = df_to_read['Max_Scr_Proc_Q_len'].median()
        self.max_scr_rep_q_len = df_to_read['Max_Scr_Rep_Q_len'].median()
        self.max_colpo_q_len = df_to_read['Max_Colpo_Q_len'].median()
        self.max_biop_proc_q_len = df_to_read['Max_Biop_Proc_Q_Len'].median()
        self.max_biop_rep_q_len = df_to_read['Max_Biop_Rep_Q_Len'].median()
        self.max_treatment_q_len = df_to_read['Max_T/t_Q_len'].median()

        self.med_gynae_res_util = df_to_read['Gynae_Res_%_util'].median()
        self.med_gynae_consul_util = df_to_read['Gynae_consul_%_util'].median()
        self.med_path_util = df_to_read['Path_%_util']
        self.med_cytotec_util = df_to_read['Cytotec_%_util']

        self.med_time_to_scr = df_to_read['Time_to_screening_results'].median()
        self.med_time_to_colpo = df_to_read['Time_to_colposcopy'].median()
        self.med_time_to_tt = df_to_read['Time_to_treatment'].median()
        self.med_tot_time_in_sys = df_to_read['Total_time_in_system'].median()



    def populate_final_summary_table(self):
        '''
        Updates the final summary table one row whenever it is called. 
        '''
        with open ('final_summary_table.csv', 'a') as f:
            writer = csv.writer(f, delimiter= ',')
            row_to_add = [parameters.experiment_no,
            self.max_scr_proc_q_len,
            self.max_scr_rep_q_len,
            self.max_colpo_q_len,
            self.max_biop_proc_q_len,
            self.max_biop_rep_q_len,
            self.max_treatment_q_len,

            self.med_gynae_res_util,
            self.med_gynae_consul_util,
            self.med_path_util,
            self.med_cytotec_util,

            self.med_time_to_scr,
            self.med_time_to_colpo,
            self.med_time_to_tt,
            self.med_tot_time_in_sys

            ]
            writer.writerow(row_to_add)
    
def clear_csv_file():
    '''f
    Erases all the contents of a csv file. Used in the refresh button of the gradio app to start fresh
    '''
    
    parameters.experiment_no = 0
    with open ('final_summary_table.csv', 'w') as f:
        pass
    open_final_table = summary_statistics()
    open_final_table.gen_final_summary_table()


def plotly_plotter():
    filepath = 'final_summary_table.csv'
    df_to_plot = pd.read_csv(filepath)
    fig = sp.make_subplots(rows = 1, cols= 3, subplot_titles= ("Max Queue Length for Different Processes", 
        'Percent Utilisation for different Professionals','Time to important events'))

    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Max_Scr_Proc_Q_len'], name = "Q len for Screen Processing"), row = 1, col = 1)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Max_Scr_Rep_Q_len'], name = "Q len for Screen Reporting"),row = 1, col = 1)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Max_Colpo_Q_len'], name = "Q len for Colposcopy"),row = 1, col = 1)        
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Max_Biop_Proc_Q_Len'], name = "Q len for Biopsy Processing"),row = 1, col = 1)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Max_Biop_Rep_Q_Len'], name = "Q len for Biopsy Reporting"),row = 1, col = 1)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Max_T/t_Q_len'], name = "Q len for Treatment"),row = 1, col = 1)
        
        
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Gynae_Res_%_util'], name = "% Utilisation for Gynae Residents"),row = 1, col = 2)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Gynae_consul_%_util'], name = "% Utilisation for Gynae Consultants"),row = 1, col = 2)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Path_%_util'], name = "% Utilisation for Pathologists"),row = 1, col = 2)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Cytotec_%_util'], name = "% Utilisation for Cytotechnicians"),row = 1, col = 2)

    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Time_to_screening_results'], name = "Time to screening results"),row = 1, col = 3)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Time_to_colposcopy'], name = "Time to Colposcopy"),row = 1, col = 3)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Time_to_treatment'], name = "Time to Treatment"),row = 1, col = 3)
    fig.add_trace(go.Scatter(x = df_to_plot['Experiment_No'], y = df_to_plot['Total_time_in_system'], name = "Total Time in the System"),row = 1, col = 3)

    return fig


def gen_kpi_table():
        #defining the KPI Results Table for one run, the export to row function in the cacx pathway class
        #  will add one row at a time
    with open('kpi_trial_results.csv', 'w') as f:
        writer = csv.writer(f, delimiter= ',')
        column_headers = [
                "Run_Number" ,
                "Max_Scr_Proc_Q_len" ,
                'Max_Scr_Rep_Q_len' ,
                'Max_Colpo_Q_len' ,
                "Max_Biop_Proc_Q_Len" ,
                "Max_Biop_Rep_Q_Len" ,
                "Max_T/t_Q_len" ,

                #resource utilisation %
                'Gynae_Res_%_util',
                'Gynae_consul_%_util',
                'Path_%_util',
                'Cytotec_%_util',

                #Time between important events
                'Time_to_screening_results',
                'Time_to_colposcopy',
                'Time_to_treatment',
                'Total_time_in_system' ]

        writer.writerow(column_headers)

open_final_table = summary_statistics()
open_final_table.gen_final_summary_table()

def main(pt_per_day, num_gynae_res, num_gynae_consul, num_cytotec, num_path):
    '''
    This function will run the simulation for different independent variables that we need.
    '''
    parameters.experiment_no += 1
    print (f'Experiment Number: {parameters.experiment_no}')
    #print('For this experiment, Pt interarrival time = 480/patients per day = 480/{pt_per_day}')
    parameters.pt_per_day = pt_per_day
    sum_stats = summary_statistics()
    
    gen_kpi_table()
    for run in range (parameters.number_of_runs):
        print(f'Run {run+1} in {parameters.number_of_runs}')
        my_sim_model = Ca_Cx_pathway(run, num_gynae_res, num_gynae_consul, num_cytotec, num_path)
        my_sim_model.run()
    sum_stats.calculate_summary_statistics()
    sum_stats.populate_final_summary_table()
    return plotly_plotter()



with gr.Blocks() as app:
    gr.HTML(
                '''
                <h1>Cervical Cancer DES App</h1>
                '''
            )
    
    with gr.Row(equal_height= True):            
                
        with gr.Column(scale = 1):
            gr.HTML(
                        '''
                        <h2>List of Assumptions</h2>
                        Discrete Event Simulation Models take in a lot of assumptions into consideration. Here is a list of assumptions made
                        while building this model
                        <ul>
                            <li> 50 Patients for cervical cancer screening arrive every day, making pt interrarrival time 480/50, but can be modified
                            <li> One working day is from 9 am to 5 pm so total 480 minutes in a day
                            <li> Procedure rooms such as Colposcopy and OT rooms are only available on Monday, Wednenday and Friday (3 days a week)
                            <li> Cytotechnicians and Pathologists process and report both Screening and Biopsy specimens
                            <li> Screening only happens through conventional pap smear
                            <li> Gynae residents and Gynae consultants essentially mean two different levels of health cadres, with Consultants being more 
                            experienced than residents. Doesn't necessarily mean professionals with specialised training in Obstetrics and Gynaecology
                            <li> In the demo app, one run is for 1000 minutes i.e. 3 days
                            <li> Medians are calculated from 100 iterations of the simulation
                        </ul>
                        '''
                    )
        with gr.Column(scale=1):
            gr.HTML('''
            <h2>AIIMS Bhopal Cervical Cancer Pathway</h2>
            ''')
            gr.Image('AIIMS_Bhopal_baseline_process_map.png')


    with gr.Row():
            gr.HTML(
                        '''
                        <h2>Modifiable Parameters</h2>
                        Limited to different HR and Procedure rooms for this implementation
                        '''
                        #num_gynae_residents, num_gynae_consultants, num_pathologists, num_cytotechnicians, num_colposcopy_room, num_ot_rooms
                    )
            pt_per_day = gr.Slider(minimum=  1, maximum = 100, label= 'Patients Visiting per day', value= 50, step = int)
            num_gynae_res = gr.Slider(minimum=  1, maximum = 10, label= 'No of Gynae Residents', value= 1, step = int)
            num_gynae_consul = gr.Slider(minimum=  1, maximum = 10, label= 'No of Gynae Consultants', value= 1, step = int)
            num_cytotec = gr.Slider(minimum=  1, maximum = 10, label= 'No of Cytotechnicians', value= 1, step = int)
            num_path = gr.Slider(minimum=  1, maximum = 10, label= 'No of Pathologists', value= 1, step = int)
    
    with gr.Row():
        btn = gr.Button(value= "Run the Simulation")
            
    with gr.Row(equal_height=True):
        output = gr.Plot(label= 'Simulation Results')
        btn.click(main, [pt_per_day, num_gynae_res, num_gynae_consul, num_cytotec, num_path], output)

    with gr.Row():
        btn_ref = gr.Button(value = "Refresh the plots")
        btn_ref.click (clear_csv_file)

app.launch(share = True)




