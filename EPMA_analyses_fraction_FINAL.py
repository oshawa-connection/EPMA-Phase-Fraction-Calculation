#!/usr/bin/env python3
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
file_count = 0
switch = 0
#Most recent imported matrix
IMPORTED_MATRIXX = np.zeros([2,2])
#Final matrix to be plotted
Final_matrix = np.zeros([1,70,1])
#Matrix to be concatonated in the third dimension
added_back = np.zeros([1,1,1])
#List of names of files for the legends.
phase_name_list=[]
#Number of bulk composition files already read
bulk_count = 0
#Size of analyses arrays
size_array = [0]
#--------------------------------------------------------------------------------------------
#-----Begin import of .tbl files.
for file in os.listdir(pwd):
    if file.endswith('_bulk.csv') and (bulk_count == 0):
    
        Bulk_composition = np.loadtxt(open(file,'rb'), dtype=float, comments='#', delimiter=',', converters=None, skiprows=1, usecols=(range(1,13)), unpack=False, ndmin=2)
        bulk_count += 1
    
    if file.endswith('_analyses.csv'):
        print ('Importing ' + file)
        IMPORTED_MATRIX = np.loadtxt(open(file,'rb'), dtype=float, comments='#', delimiter=',', converters=None, skiprows=1, usecols=(range(1,13)), unpack=False, ndmin=2)
        print ('Successfully imported ' + file)

        IMPORTED_MATRIXX = np.zeros([(IMPORTED_MATRIX.shape[0]),(IMPORTED_MATRIX.shape[1])])
        IMPORTED_MATRIXX = IMPORTED_MATRIX

        phase_name_list.append(file)
        
        file_count += 1
        switch = 1
	#Set the final matrix first z value equal to the imported file matrix.

    if (switch == 1) and (file_count == 1):
        #What dimensions
        rows = IMPORTED_MATRIXX.shape[0]
        cols = IMPORTED_MATRIXX.shape[1]
        Final_matrix = np.zeros([rows,cols,1])
        Final_matrix[:,:,0] = IMPORTED_MATRIXX
        size_array[0] = rows
        #Set Switch back to 0
        switch = 0

    #If the newly imported matrix has more rows than the previous rows, we must then extend the previous matrix with nan's
    elif (switch ==1) and ((IMPORTED_MATRIXX.shape[0])>(Final_matrix.shape[0])):
    
        new_matrix = np.zeros([IMPORTED_MATRIXX.shape[0],IMPORTED_MATRIXX.shape[1],file_count])
        new_matrix.fill(np.nan)
        for i in range(file_count):
            new_matrix = inner_product_two(new_matrix, Final_matrix[:,:,i-1],i)
                        
        new_matrix[:,:,(file_count-1)] = IMPORTED_MATRIXX

        Final_matrix = np.zeros([IMPORTED_MATRIXX.shape[0],IMPORTED_MATRIXX.shape[1],file_count])
        Final_matrix = new_matrix
        
        #Extend the array which will be used for the loop later.
        size_array.extend([IMPORTED_MATRIX.shape[0]])
        
        switch = 0
    #If the newly imported matrix is smaller or the same size as the storage matrix, then just add it to the newest dimension.
    #As normal.
    #If the file is not a .csv then skip it.
    elif switch==1:

        added_back = np.zeros([Final_matrix.shape[0],Final_matrix.shape[1],1])
        added_back.fill(np.nan)

        Final_matrix = np.append(Final_matrix,added_back,2)
        Final_matrix = inner_product_two(Final_matrix,IMPORTED_MATRIXX,(file_count))
        #Set Switch back to 0
        switch = 0
        size_array.extend([IMPORTED_MATRIX.shape[0]])

#Make sure the user actually provided an input file.
if file_count == 0:
    print ('You need to provide an _analyses.csv input file.')
    sys.exit()
else:
    print ('Total number of files read and successfully concatonated = %i'%file_count) 
    Analyses_matrix = Final_matrix 
    
#----------------------------------------------------------------------------------------------------------------------------
#Number of phases
m = file_count 

#Number of possible combinations of analyses.
number_analyses_combos = 1
for i in range(0,m):
    number_analyses_combos = number_analyses_combos * size_array[i]

#Number of equations
#Number of possible combinations of analyses
n = (number_analyses_combos * Analyses_matrix.shape[1]) + 1

#Number of phases to consider.
m = Analyses_matrix.shape[2]

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
x = np.linalg.lstsq(a,b)
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