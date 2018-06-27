#!/usr/local/bin/python3.6
# -*- coding: utf-8 -*-
"""Created on Fri Feb 23 18:35:38 2018
Fitting routine for calculating the weight fraction of phases from bulk rock 
chemistry and EPMA analysis of solid phase composition.

Problems with the calculation are that if phases have very similar phase
composition, then they cannot easily be distinguished.
"""
#--------------------------------------------------------------------------------------------
#-----Import all modules used
import os
import numpy as np
import pandas as pd
from itertools import product
#------------------------------------------------------------------------------
#----Load in all analyses
#-----Initialise
#List all files within the current directory.
pwd = os.getcwd()
#List of names of files for the legends.
phase_name_list=[]
#Number of bulk composition files already read
bulk_count = 0
#Size of analyses arrays
size_array = [0]
found = False
#--------------------------------------------------------------------------------------------
#-----Import .xls files.

analyses = pd.read_excel('analyses.xlsx')
analyses = pd.read_excel('analyses.xlsx')
analyses.drop(['comment'],axis=1,inplace=True)


grouped_analyses = analyses.groupby('phase')


#----------------------------------------------------------------------------------------------------------------------------
#Number of phases
n = len(grouped_analyses)

#Number of possible combinations of analyses.
number_analyses_combos = 1
for i in range(0,n):
    number_analyses_combos = number_analyses_combos * size_array[i]


#Number of possible combinations of analyses
m = (number_analyses_combos * Analyses_matrix.shape[1]) + 1


analyses.drop(['comment'],axis=1,inplace=True)

sizes_array = []
phase_list = []


for phase, phase_df in grouped_analyses:
    if phase == 'bulk':
        Bulk_composition = np.array(phase_df.drop(['phase'],axis=1))
    else:
       # sizes_array.append(len(phase_df))
        phase_list.append(phase)
        
    
    
    

    



#Equation
#a . x  = b
#------------------------------------------------------------------------------
#----Initialise Matrices
#The initial bulk composition of the experiment
b = np.zeros([m,1]) #b = n * 1

#The coefficient matrix of fraction components for each analysis of each phase
a = np.zeros([n,m]) #a = n * m
#The vector of weight percentages of each phase
x = np.zeros([n,1]) #x = m * 1

sizes_array = []

for i in range(0,m):
    sizes_array.append(range(size_array[i]))

#Get all possible combinations of the phase analyses, encoded as the analysis number
#And stored in an array.
#Also take LSD and move out to the jungle somewhere     
idm = list(product(*sizes_array))
idm = np.asarray(idm)   
row_counter = 0

#Fill Coefficient matrix.
for component in range(0,Analyses_matrix.shape[1]):
    for i in range(0,idm.shape[0]):
        for j in range(0,idm.shape[1]):
            a[row_counter,j] = Analyses_matrix[idm[i,j],component,j]
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
if any(i<0 for i in x) == True:
    print('')
    print('** Solved for negative values, switching to different solver ** ')
    print('** Solver may give weird values. **')
    from scipy.optimize import nnls
    b = b[:,0]
    x = nnls(a, b)
    residuals = x[1]
    print('Solved. Residual is %f' %(residuals))
    x = x[0]
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
