# EPMA Phase Fraction Calculation
Calculate the weight percentage of each phase in a system from bulk system composition and phase analyses.

If the bulk composition of the system (rock, magma) is known, and point analyses of all phases has been conducted using EMPA/EPMA or another technique,  the weight fraction that each phase contributes to the bulk system can be calculated using a least squares regression fitting technique. 

This avoids using so called NORM techniques which are complicated, make several simplifying assumptions about the partitioning of elements into certain phases, and are limited in which phases are reported to be present.

To use this file, you need python 2/3 and numpy installed. Provide one file with something_bulk.csv as the filename to represent the bulk system composition. Provide as many other files with the title phase1_analyses.csv as there are phases in the system, examples included. 

Then run the script using python. Output is as text in the terminal.
