# -*- coding: utf-8 -*-
"""
Created on 18 October 2024 14:57 (HSC DPhil students Room, Oxford)

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


#Non Modifiable variables
class parameters (object):
    '''
    This class contains all the constant non-modifiable parameters that will go into the model
    These mostly include service times and the time for which the simulation is supposed to run and halt etc
    '''

    pt_per_day = None #Number of patients that visit the Gynae OPD every day (derived from AIIMS Bhopal Annual Report)(This could also be capped for a day if there are limited spots)
    pt_interarrival_time = None #Number of minutes in a working day / number of total patients expected during the day
    #modifiable factors, will be defined again in the relevant class
    #resources 
    #Staff
    gynae_resident = None #Gynaecological residents that perform history and examination and pap smear collection in the routine OPDs in AIIMS Bhopal
    gynae_consulant = None #Gynaecological consultants that perform procedures such as LEEP (LEETZ) or hysterectomy etc
    pathologists = None #No of pathologists that interpret the pathology findings
    cytotechnicians = None #No of cytotechnicians that process the sample generated


    #Stuff
    num_pap_kits = None #number of kits that the hospital has to perform a pap smear (ayre's spatula, glass slide, preservative and box)
    num_pathology_consumables = None #consumables required for processing the pathological specimen
    num_colposcopy_consumables = None #consumables required for conducting colposcopy
    num_thermal_consumables = None #consumables required for LEEP
    num_ot_consumables = None #consumables required for hysterectomy


    #rooms (scheduled resource)
    num_colposcopy_rooms = None #number of colposcopy rooms
    num_ot_rooms = None #number of OT rooms

    #service times
    history_exam_time = None #time taken to complete history and examination per patient (Imp thing to remember here would be that this might change as the system adapts
                            #to an excess load)
    path_processing_time = None #time it takes from the sample is generated to the sample is prepared by cytotechnicians and ready for interpretation
    path_reporting_time = None # time taken by pathologists to report the results of a processed sample
    colposcopy_time = None #time taken to perform 1 colposcopy
    thermal_time = None #time taken for 1 Loop Electrosurgical Excision procedure
    hysterectomy_time = None #time taken for 1 hysterectomy

    #Epidemiological parameters
    screen_positivity_rate = None #% of positive samples (True positive + false positive / total samples)
    biopsy_rate = None # % of all colposcopies that undergo a biopsy
    biopsy_cin_rate = None # % of biopsies that are CIN
    biopsy_cacx_rate = None # % of biopsies that are CaCx
    follow_up_rate = None # % of women who follow up after a positive screen result (My own meta analysis + local data)

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
        self.pt_id = pt_id
        #declaring the variables to be recorded
        self.time_at_entered = None #time when the patient entered into the OPD room
        self.time_at_screen_result = None #time when the patient first received the screening result
        self.time_at_colposcopy = None #time when the patient attended the colposcopy clinic
        self.time_at_treatment = None #time when patient got the treatment, either admission or surgery or LEEP or thermal/cryo
        self.history_examination_service_time = None
        self.colposcopy_service_time = None
        self.treatment_service_time = None

class screen_sample(object):
    '''
    This class creates an instance of a screening sample object that undergoes processing and reporting
    '''
    def __init__(self, ss_id):
        self.screen_sample_id = ss_id #sample id is the same as the patient id, when creating an instance of the sample id, make sure to enter patient id 

        #declaring variables that will be recorded later on
        self.screen_sample_processing_time = None
        self.screen_sample_reporting_time = None


class biopsy_sample(object):
    '''
    This class generates an instance of a biopsy sample object that undergoes processing and reporting
    '''
    def __init__(self, bs_id):
        self.biopsy_sample_id = bs_id #sample id is the same as the patient id, when creating an instance of the biopsy sample id, 
                                        #make sure to enter patient id 

        self.biopsy_sample_processing_time = None
        self.biopsy_sample_reporting_time = None


class Ca_Cx_pathway (object):
    '''
    This is the fake hospital. Defines all the processes that the patients will go through. Will record statistics for 1 simulation that will be later analyzed and clubbed with 
    results from 100 simulations
    '''
    def __init__(self, num_gynae_residents, num_gynae_consultants, num_pathologists, num_cytotechnicians, num_colposcopy_room, num_ot_rooms):
        self.env = simpy.Environment()

        #declaring number of modifiable resource capacity, non modifiable resources to be imported from the parameters class
        self.num_gynae_residents = num_gynae_residents
        self.num_gynae_consultants = num_gynae_consultants
        self.num_pathologists = num_pathologists
        self.num_cytotechnicians = num_cytotechnicians

        self.colposcopy_schedule = [] #list of integers form 0-6 for each day of the week that resource is available
        self.ot_schedule = []  #list of integers from 0-6 for each day of the week that resource is available
        
        self.pt_counter = None #acts as the UHID of the patient
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
            "Time_at_screening_result":[],
            "Time_at_colposcopy" : [],
            "Time_at_treatment" : [],
            "History and Examination time": [], #also recording service times as they will ultimately be added up to calculate resource utilisation percentage
            "screen_processing_time":[],
            "Screen_reporting_time":[],
            "Colposcopy_time":[],
            "Treatment_time":[],
        })

        #declaring system KPIs to be measured at the run level.
        #Queue lengths for different processes
        self.max_q_len_screen_processing = None
        self.max_q_len_screen_reporting = None
        self.max_q_len_colposcopy = None
        self.max_q_len_colposcopy_processing = None
        self.max_q_len_colposcopy_reporting = None
        self.max_q_len_treatment = None

        #Resource utilization percentages
        self.gynae_residents_utilisation = None
        self.gynae_consultants_utlisation = None
        self.cytotechnician_utilisation = None
        self.pathologist_utilisation = None
    
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
            if self.is_within_working_hours:
        #if time of day is appropriate then generate the patient
                self.patient = ca_cx_patient(self.pt_counter)
                self.pt_counter += 1
        #record necessary timepoints
                self.patient.time_entered = self.env.now
        #patient moves to the OPD
                self.env.process(self.history_examination(self.patient))
        #time for next patient arrival
                wait_time_for_next_pt = random.expovariate(1/parameters.pt_interarrival_time)
                yield self.env.timeout(wait_time_for_next_pt)

    def history_examination(self, patient):
        '''
        Patient undergoes history and examination and in the process also generates the screening sample
        '''
        #request for a resident and consumables for sample collection and wait for them to be available
        with self.gynae_residents.request() as gynae_res, self.pap_kit.request() as pap, self.screening_consumables as scr_consum :
            yield gynae_res and pap and scr_consum
        
        #patient undergoes history, examination and sample collection
            history_examination_time = random.triangular(parameters.history_exam_time/2, parameters.history_exam_time, parameters.history_exam_time *2 )
            yield self.env.timeout(history_examination_time)

            #generate a screening sample
            self.screening_sample = screen_sample(patient.pt_id) #screen sample id is the same as the patient id
            #sample goes on for processing
            self.env.process(self.screen_sample_processing(self.screening_sample))


    def screen_sample_processing(self, screening_sample):
        '''
        Sample undergoes processing
        '''
        #request resources and wait for them to be available
        with self.cytotechnician.request() as cytotec, self.pathology_consumables.request() as scr_proc_consum:
            yield cytotec and scr_proc_consum
        
        #sample undergoes processing
            screen_sample_processing_time = random.triangular(parameters.path_processing_time/2, parameters.path_processing_time, parameters.path_processing_time *2)
            yield self.env.timeout(screen_sample_processing_time)
        #sample goes for reporting
            self.env.process(self.screen_sample_reporting(self.screening_sample))

    def screen_sample_reporting(self, screening_sample):
        '''
        Processed sample is interpreted and reported by pathologist
        '''
        #request for a pathologist and wait until 
        with self.pathologist.request() as path:
            yield path

        #record the current time as an important milestone
            self.patient.time_at_screen_result = self.env.now
        
        #sample undergoes reporting
            screen_sample_reporting_time = random.triangular(parameters.path_reporting_time/2, parameters.path_reporting_time, parameters.path_reporting_time * 2)
            yield self.env.timeout(screen_sample_reporting_time)

        #if sample is positive, move on to follow up, otherwise terminate
            if random.random() < parameters.screen_positivity_rate:
                self.env.process(self.call_for_follow_up())

    def call_for_follow_up (self, screen_sample):
        '''
        Gynaecology residents 
        '''
        #request a gynae_res (later on could modify to include a receptionist or another health cadre)
        with self.gynae_residents.request() as gynae_res:
            yield gynae_res
        
        # whether the patient returns or not
            if random.random < parameters.follow_up_rate:
            #patient goes on for colposcopy
                self.env.process(self.colposcopy(self.patient))
        #instantaneous process so no timeout really and also not a service

    def colposcopy(self):
        '''
        Patient that was generated undergoes colposcopy
        '''
        #requests for a consultant, consumables and a room
        with self.gynae_consultants.request() as gynae_consul, self.colposcopy_consumables.request() as gynae_consumables, self.colposcopy_room.request() as colpo_room:
            yield gynae_consul and gynae_consumables and colpo_room

        #Record time at colposcopy
            self.patient.time_at_colposcopy = self.env.now

        #check if biopsy is performed or not
            if random.random() < parameters.biopsy_rate:
                pt_biopsy_sample = biopsy_sample(self.patient.pt_id) #generate a biosy sample that will go for processing  
                #biosy sample goes for processing
                self.env.process(self.biopsy_sample_processing(pt_biopsy_sample))
        


    def biopsy_sample_processing(self, biopsy_sample):
        '''
        Biopsy sample if prepared undergoes processing
        '''
        #requests a cytotechnicians and consumables
        with self.cytotechnician.request() as cytotec, self.pathology_consumables.request() as path_consum:
            yield cytotec and path_consum

            #biopsy sample undergoes processing
            biopsy_sample_processing_time = random.triangular(parameters.path_processing_time/2, parameters.path_processing_time, parameters.path_processing_time * 2)    
            yield self.env.timeout(biopsy_sample_processing_time)

            #biopsy sample goes for reporting
            self.env.process(self.biopsy_sample_reporting(biopsy_sample))


    def biopsy_sample_reporting(self, biopsy_sample):
        '''
        Biopsy sample if taken undergoes reporting after processing
        '''
        #requests a pathologist
        with self.pathologist.request() as path:
            yield path

            #biopsy sample undergoes reporting
            biopsy_sample_reporting_time = random.triangular(parameters.path_reporting_time/2, parameters.path_reporting_time, parameters.path_reporting_time *2)
            yield self.env.timeout(biopsy_sample_reporting_time)

            #depending on the diagnosis, patient either goes for thermal ablation or hysterectomy (currently, only making 2 options available, have the option of adding more on later)
            biopsy_result = random.random()
            if biopsy_result < parameters.biopsy_cin_rate:
                self.env.process(self.thermal_ablation(self.patient)) #diagnosed with CIN
            
            elif parameters.biopsy_cin_rate < biopsy_result < parameters.biopsy_cacx_rate:
                self.env.process(self.hysterectomy(self.patient)) #diagnosed with cervical cancer
            
            else:
                self.patient.time_at_treatment = self.env.now #patient exits the system

    def thermal_ablation(self, patient):
        '''
        If indicated, pt undergoes thermal ablation
        '''
        #requests resources required for thermal ablation
        with self.gynae_consultants.request() as gynae_consul, self.thermal_consumables.request() as thermal_consum, self.colposcopy_room as colpo_room:
            yield gynae_consul and thermal_consum and colpo_room

        #patient undergoes thermal ablation
            thermal_ablation_time = random.triangular(parameters.thermal_time/2, parameters.thermal_time, parameters.thermal_time *2)
            yield self.env.timeout(thermal_ablation_time)

        #patient exits the system
            self.patient.time_at_treatment = self.env.now

    def leep (self,patient):
        '''
        if indicated, patient undergoes LEEP
        '''
        #Not being implemented in this first version of the model
        pass 
    def hysterectomy (self,patient):
        '''
        if indicated, patient undergoes hysterectomy
        '''
        #request for a ot room and other equipment
        with self.gynae_consultants.request() as gynae_consul, self.ot_consumables.request() as ot_consum, self.ot_room as ot_room:
            yield gynae_consul and ot_consum and ot_room

            #patient undergoes surgery
            hysterectomy_time = random.triangular(parameters.hysterectomy_time/2, parameters.hysterectomy_time, parameters.hysterectomy_time *2)
            yield self.env.timeout(hysterectomy_time)

            #patient exits the system
            






class summary_statistics(object):
    '''
    This class will define methods that will calculate aggregate statistics from 100 simulations and append the results onto a new spreadsheet which will be used to append results
    from 100 simulations for different number of independent variables (such as patients)
    '''

class Plotter (object):
    '''
    This class will define methods that will plot data collated into various spreadsheets in the summary statistics class and generate nice plots. Which will either be saved 
    on the local machine using matplotlib or be displayed on a web browser using plotly
    '''

class web_app (object):
    '''
    This class will define a gradio application that will contain the diagram of the process map and ability for the stakeholders to modify the different parameters
    and conduct experiments themselves.
    '''

# def main():
'''
This function will run the simulation for different independent variables that we need.
'''

#Goal will be to complete one class at a time, not necessarily in one sitting, but each commit will have one class or so.
#This is the skeleton on and the initial commit
