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


#Non Modifiable variables
class constants (object):
    '''
    This class contains all the constant non-modifiable parameters that will go into the model
    These mostly include service times and the time for which the simulation is supposed to run and halt etc
    '''

class patients (object):
    '''
    This class creates patients and declares their individual parameters that explains how they spent their time at the hospital
    These individual parameters will then be combined with others in the simulation to get overall estimates
    '''

class ED_sim (object):
    '''
    This is the fake hospital. Defines all the processes that the patients will go through. Will record statistics for 1 simulation that will be later analyzed and clubbed with 
    results from 100 simulations
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
