#!/usr/local/bin/python3.6
# -*- coding: utf-8 -*-
"""Created on Fri Feb 23 18:35:38 2018
Fitting routine for calculating the weight fraction of phases from bulk rock 
chemistry and EPMA analysis of solid phase composition.

Problems with the calculation are that elements that have a high
weight are favoured over those with small weights.
If the calculation was done with molar fraction,
then crazy components like Si4O8 and Al4O6 instead of SiO2 and Al2O3, could 
be used to bring them into similar ranges as each other again.

"""
#--------------------------------------------------------------------------------------------
#-----Import all modules used
import sys
import os
import numpy as np
import pandas as pd
from itertools import product
#------------------------------------------------------------------------------
#----Load in all analyses
#-----Initialise
#List all files within the current directory.
pwd = os.getcwd()
#Size of analyses arrays
size_array = [0]
found = False
#--------------------------------------------------------------------------------------------
#-----Begin import of .tbl files.
for filename in os.listdir(pwd):
    if filename.lower() == ('analyses.xlsx'):
        analyses = pd.read_excel(filename)
        found = True
    elif filename.lower() == 'analyses.csv':
        analyses = pd.read_csv(filename)
        found = True
    
if found != True:
    print ('You need to provide an analyses.csv or analyses.xlsx input file.')
    sys.exit()

#----------------------------------------------------------------------------------------------------------------------------
#Number of phases



#m = len(phase_list)
#
##Number of possible combinations of analyses.
#number_analyses_combos = 1
#for i in range(0,m):
#    number_analyses_combos = number_analyses_combos * size_array[i]
#
##Number of equations
##Number of possible combinations of analyses
#n = (number_analyses_combos * Analyses_matrix.shape[1]) + 1


analyses.drop(['comment'],axis=1,inplace=True)

grouped_analyses = analyses.groupby('phase')
phase_list = []
analyses_matrix_list = []
bulk_switch = False
max_size = 0

for phase, phase_df in grouped_analyses:
    if phase.lower() == 'bulk':
        if bulk_switch == False:
        
            Bulk_composition = np.array(phase_df.drop(['phase'],axis=1))
            bulk_switch = True
    else:
        
        phase_list.append(phase)
        
        analyses_matrix_list.append(np.array(phase_df.drop(['phase'],axis=1)))
                
        if phase_df.drop(['phase'],axis=1).shape[0] > max_size:
            max_size = phase_df.drop(['phase'],axis=1).shape[0] 
        
    
#Make sure bulk entry was added.    
if bulk_switch == False:
    print ('Make sure you provided a phase called bulk in the .xls file.')
    sys.exit()

#Number of components to consider
n_components = analyses_matrix_list[0].shape[1]
iteration = 1
for phase_array in analyses_matrix_list:    
    vstack_size = [(max_size - (phase_array.shape[0])),n_components]    
    vstack_array = np.zeros(vstack_size)
    vstack_array.fill(np.nan)
    #Vstack to extend with NaN
    dstack_array = np.vstack((phase_array,vstack_array))    
    #Dstack to add to 3rd dimension
    if iteration == 1:
        final_matrix = dstack_array
        iteration += 1
    else:        
        final_matrix =  np.dstack((final_matrix,dstack_array))
    
    
sys.exit()
#Equation
#a . x  = b
#------------------------------------------------------------------------------
#----Initialise Matrices
#The initial bulk composition of the experiment
b = np.zeros([n,1]) #b = n * 1

#The coefficient matrix of fraction components for each analysis of each phase
a = np.zeros([n,m]) #a = n * m
#The vector of weight percentages of each phase
x = np.zeros([m,1]) #x = m * 1

sizes_array = []

for i in range(0,m):
    sizes_array.append(range(size_array[i]))

#Get all possible combinations of the phase analyses, encoded as the analysis number
#And stored in an array.
#ACM = analyses combinations matrix.
acm = list(product(*sizes_array))
acm = np.asarray(acm)   
row_counter = 0

#Fill Coefficient matrix.
for component in range(0,Analyses_matrix.shape[1]):
    for i in range(0,acm.shape[0]):
        for j in range(0,acm.shape[1]):
            a[row_counter,j] = Analyses_matrix[acm[i,j],component,j]
        row_counter +=1

#Sum of all components is  = 1. All coefficients are 1 also.                       
a[row_counter,:] = 1

#Fill the RHS vector with bulk composition analyses.
counter2 = 0
for j in range(0,Analyses_matrix.shape[1]):
    for i in range(0,number_analyses_combos):
        b[counter2,0] = Bulk_composition[0,j]
        counter2 += 1 

#Sum of all components is  = 1
b[n-1,0] = 1
#Solve by using least squares regression, to find the best fit
x,residuals,rank,s = np.linalg.lstsq(a,b)
# =============================================================================
# Calculate and print residuals.
# =============================================================================
r2 = 1 - residuals / (b.size * b.var())
print ('Success. R^2 is: %f'%(r2))

x = x*100

total = np.sum(x)

#Print output
print ('')
print ('Total weight fraction of phases are %f %%'%total)
for i in range(0,len(phase_name_list)):
    phase_name_list[i] = (phase_name_list[i])[:-13] 
    print ('Weight fraction of %s is %f %%' % (phase_name_list[i],x[i]))

#Normalise and reprint
x = x / total * 100

print ('')
print ('Total, normalised weight fraction of phases are 100 %')
for i in range(0,len(phase_name_list)):  
    print ('Weight fraction of %s is %f %%' % (phase_name_list[i],x[i]))   
