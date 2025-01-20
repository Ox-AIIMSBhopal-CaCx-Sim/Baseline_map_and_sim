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
 
    #modifiable factors, will be defined again in the relevant class
    #resources 
    #Staff
    gynae_resident = None #Gynaecological residents that perform history and examination and pap smear collection in the routine OPDs in AIIMS Bhopal
    gynae_consulant = None #Gynaecological consultants that perform procedures such as LEEP (LEETZ) or hysterectomy etc
    pathologists = None #No of pathologists that interpret the pathology findings
    cytotechnicians = None #No of cytotechnicians that process the sample generated


    #Stuff
    num_pap_kits = None #number of kits that the hospital has to perform a pap smear (ayre's spatula, glass slide, preservative and box)
    num_screening_consumables = None #consumables required for screening
    num_pathology_consumables = None #consumables required for processing the pathological specimen
    num_colposcopy_consumables = None #consumables required for conducting colposcopy
    num_leep_consumables = None #consumables required for LEEP
    num_ot_consumables = None #consumables required for hysterectomy


    #rooms (scheduled resource)
    num_colposcopy_rooms = None #number of colposcopy rooms
    num_ot_rooms = None #number of OT rooms

    #service times
    history_exam_time = None #time taken to complete history and examination per patient (Imp thing to remember here would be that this might change as the system adapts
                            #to an excess load)
    sample_processing_time = None #time it takes from the sample is generated to the sample is prepared by cytotechnicians and ready for interpretation
    sample_reporting_time = None # time taken by pathologists to report the results of a processed sample
    colposcopy_time = None #time taken to perform 1 colposcopy
    leep_time = None #time taken for 1 Loop Electrosurgical Excision procedure
    hysterectomy_time = None #time taken for 1 hysterectomy


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
        self.time_entered = None #time when the patient entered into the OPD room
        self.time_screen_result = None #time when the patient first received the screening result
        self.time_colposcopy = None #time when the patient attended the colposcopy clinic
        self.time_treatment = None #time when patient got the treatment, either admission or surgery or LEEP or thermal/cryo


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
        
        #declaring resources
        #staff
        self.gynae_residents = simpy.Resource(self.env, capacity=num_gynae_residents)
        self.gynae_consultants = simpy.Resource(self.env, capacity=num_gynae_consultants)
        self.pathologist = simpy.Resource(self.env, capacity=num_pathologists)
        self.cytotechnician = simpy.Resource(self.env, capacity=num_cytotechnicians)

        #stuff
        self.pap_kit = simpy.Resource(self.env, capacity=parameters.num_pap_kits)
        self.screening_consumables = simpy.Resource(self.env, capacity=parameters.num_screening_consumables)
        self.pathology_consumables = simpy.Resource(self.env, capacity=parameters.num_pathology_consumables)
        self.colposcopy_consumables = simpy.Resource(self.env, capacity=parameters.num_colposcopy_consumables)
        self.leep_consumables = simpy.Resource(self.env, capacity=parameters.num_leep_consumables)
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
            "Time_at_treatment" : []
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
    
    def gen_patient_arrival(self):
        '''
        Generates a fictional patient according to a distribution, they undergo and OPD, this generates a sample which undergoes processing, after results are
        conveyed, if positive, patient only then moves on to the next step i.e. colposcopy.
        '''


    def gen_screen_sample(self):
            '''
            Generates a screening sample, in this case a pap smear. This undergoes processing and reporting, after results are conveyed, 
            '''

    def screen_sample_processing(self):
            '''
            Sample undergoes processing
            '''

    def screen_sample_reporting(self):
        '''
        Processed sample is interpreted and reported by pathologist
        '''

    def colposcopy(self, patient):
        '''
        Patient that was generated undergoes colposcopy
        '''

    def biopsy_sample_processing(self):
        '''
        Biopsy sample if prepared undergoes processing
        '''
    
    def biopsy_sample_reporting(self):
        '''
        Biopsy sample if taken undergoes reporting after processing
        '''

    def thermal_ablation(self, patient):
        '''
        If indicated, pt undergoes thermal ablation
        '''

    def leep (self,patient):
        '''
        if indicated, patient undergoes LEEP
        '''

    def hysterectomy (self,patient):
        '''
        if indicated, patient undergoes hysterectomy
        '''







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
