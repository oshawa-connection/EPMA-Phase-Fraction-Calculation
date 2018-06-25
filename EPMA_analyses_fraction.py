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
#----Define all functions.
#Matrix inner product
#Some weird quirk with numpy means that I can't do this the normal way.
def inner_product(big_matrix, small_matrix):
        
	for i in range(small_matrix.shape[0]):
		for j in range(small_matrix.shape[1]):
			big_matrix[i,j] = small_matrix[i,j] #Don't know why this doesn't work normally but oh well, this is going to be slow :(
	return big_matrix

def inner_product_two(big_matrix, small_matrix,file_count):       
	for i in range(small_matrix.shape[0]):
		for j in range(small_matrix.shape[1]):
			big_matrix[i,j,file_count-1] = small_matrix[i,j] #Don't know why this doesn't work normally but oh well, this is going to be slow :(
	return big_matrix 
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
m = len(phase_list)

#Number of possible combinations of analyses.
number_analyses_combos = 1
for i in range(0,m):
    number_analyses_combos = number_analyses_combos * size_array[i]

#Number of equations
#Number of possible combinations of analyses
n = (number_analyses_combos * Analyses_matrix.shape[1]) + 1


analyses.drop(['comment'],axis=1,inplace=True)

grouped_analyses = analyses.groupby('phase')
phase_list = []


for phase, phase_df in grouped_analyses:
    if phase == 'bulk':
        Bulk_composition = np.array(phase_df.drop['phase'])
    else:
    
        phase_list.append(phase)
    
    
    print(phase_df)
    

    



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
