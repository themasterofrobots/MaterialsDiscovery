from __future__ import division
import math
import json
import os.path

#USER INPUT
#num_k_points = 50 # Zn0 OUTCAR
#num_bands = 36
#N_CB = 2
#N_VB = 4
#energy_index_start = 8669 # Lowest energy VB line number (first line of file is line 0)
#k_point_index_start = 10716 # First k-point line number (first line of file is line 0)
#filename = 'OUTCAR.wz-ZnO'

#num_k_points = 131
#num_bands = 64
#N_CB = 2
#N_VB = 4
#energy_index_start = 12011 # Lowest energy VB line number (first line of file is line 0)
#k_point_index_start = 20966 # First k-point line number (first line of file is line 0)
#filename = 'OUTCAR.wz-GaN'

#num_k_points = 145
#num_bands = 20
#N_CB = 1
#N_VB = 2
#energy_index_start = 23932 # Lowest energy VB line number (first line of file is line 0)
#k_point_index_start = 27467 # First k-point line number (first line of file is line 0)
#filename = 'OUTCAR.rs-MgO'

#num_k_points = 50
#num_bands = 16
#N_CB = 2
#N_VB = 4
#energy_index_start = 9041 # Lowest energy VB line number (first line of file is line 0)
#k_point_index_start = 10100 # First k-point line number (first line of file is line 0)
#filename = 'OUTCAR.wz-AlN'

#num_k_points = 72
#num_bands = 32
#N_CB = 1
#N_VB = 2
#energy_index_start = 13835 # Lowest energy VB line number (first line of file is line 0)
#k_point_index_start = 16469 # First k-point line number (first line of file is line 0)
#filename = 'OUTCAR.zb-GaN'

num_k_points = 145
num_bands = 18
N_CB = 1
N_VB = 2
energy_index_start = 23942 # Lowest energy VB line number (first line of file is line 0)
k_point_index_start = 27185 # First k-point line number (first line of file is line 0)
filename = 'OUTCAR.rs-CdO'

# PARSE FILE
lines_between_k_points = num_bands + 4
k_point_index_end = k_point_index_start + num_k_points
k_point_weights = []
f = open(filename, 'r')
lines=f.readlines()
k_point_data = lines[k_point_index_start:k_point_index_end]
for j in range(k_point_index_end-k_point_index_start):
	k_point_weights.append(int(k_point_data[j][34:37])) # Creates list of all k-point weights

k_point_energy = []
weighted_energy = []
vbm = []
energy_index_end = energy_index_start + N_VB + N_CB
for n in range(num_k_points): # gets QP-energies from OUTCAR file
	x = []
	y = []
	energy_data = lines[energy_index_start+(n*lines_between_k_points):energy_index_end+(n*lines_between_k_points)]
	for k in range(energy_index_end-energy_index_start):
		x.append(float(energy_data[k][25:33]))
		y.append(float(energy_data[k][25:33])*k_point_weights[n])
		if k == N_VB-1:
			vbm.append(float(energy_data[k][25:33]))
	k_point_energy.append(x) # list of all energies
	weighted_energy.append(y)

k_point_sum_matrix = []
for u in range(num_k_points): # I think something is wrong with the math here
	sum_VB = (math.fsum(weighted_energy[u][0:N_VB]))/N_VB
	sum_CB = (math.fsum(weighted_energy[u][N_VB:(N_VB+N_CB)]))/N_CB
	k_point_sum_matrix.append(sum_VB + sum_CB)

k_point_sum = math.fsum(k_point_sum_matrix)
E_BP = (k_point_sum)/(2*math.fsum(k_point_weights)) - max(vbm)
print E_BP 

