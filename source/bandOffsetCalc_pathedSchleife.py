from __future__ import division
import math
import json
import os.path

num_k_points = 188
num_bands = 36
N_CB = 1
N_VB = 2
energy_index_start = 1943 # Lowest energy VB line number (first line of file is line 0)
filename = 'OUTCAR.bands'

# PARSE FILE
lines_between_k_points = num_bands + 3
f = open(filename, 'r')
lines=f.readlines()

k_point_energy = []
weighted_energy = []
vbm = []
energy_index_end = energy_index_start + N_VB + N_CB
for n in range(num_k_points): # gets QP-energies from OUTCAR file
	x = []
	energy_data = lines[energy_index_start+(n*lines_between_k_points):energy_index_end+(n*lines_between_k_points)]
	for k in range(energy_index_end-energy_index_start):
		x.append(float(energy_data[k][12:20]))
		if k == N_VB-1:
			vbm.append(float(energy_data[k][12:20]))
	k_point_energy.append(x) # list of all energies

k_point_sum_matrix = []
for u in range(num_k_points): # I think something is wrong with the math here
	sum_VB = (math.fsum(k_point_energy[u][0:N_VB]))/N_VB
	sum_CB = (math.fsum(k_point_energy[u][N_VB:(N_VB+N_CB)]))/N_CB
	k_point_sum_matrix.append(sum_VB + sum_CB)

k_point_sum = math.fsum(k_point_sum_matrix)
E_BP = (k_point_sum)/(2*num_k_points) - max(vbm)
print E_BP 

