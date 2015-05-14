from __future__ import division
import math
import json
import os.path
from pymatgen.electronic_structure.bandstructure import BandStructure

from multiprocessing import Pool
NUM_PROCESSES = 2

PRINT_DEBUG = True

ACTIVE_LAYER_SUBDIR = 'activeLayer'
OUTER_LAYER_SUBDIR = 'blockingLayers'

N_CB = 1
N_VB = 2

def add_ebps(dirname):
	#p = Pool(NUM_PROCESSES)
	#p.map(add_ebp, [dirname + '/' + dir_entry for dir_entry in os.listdir(dirname)])
	for dir_entry in os.listdir(dirname):
		add_ebp(dirname + '/' + dir_entry)

def add_ebp(dir_entry):
	json_data = {}
	with open(dir_entry, 'r') as thefile:
		json_data = json.load(thefile)
		if 'ebp' in json_data:
			return
		bs = BandStructure.from_dict(json_data['bandstructure'])
		for i in range(bs.nb_bands):
			if min(bs.bands[1][i]) > bs.efermi:
				cbbottom = i
				vbtop = i-1
				break

		vb_en = 0.0 
		cb_en = 0.0 
		for i in range(N_VB):
			vb_en += math.fsum(bs.bands[1][vbtop-i])/N_VB
		for i in range(N_CB):
			cb_en += math.fsum(bs.bands[1][cbbottom+i])/N_CB 

		total_en = (vb_en + cb_en)/(2*len(bs.kpoints))
		if not(bs.is_metal()):
			total_en -= (bs.get_vbm())['energy']
			json_data['ebp'] = total_en
		else:
			print 'Warning, material is a metal! No VBM offset applied, no file written'
		if PRINT_DEBUG:
			print 'The branch point energy of ' + json_data['name'] + ' is ' + str(total_en) + ' eV'
	with open(dir_entry, 'w') as thefile:
		json.dump(json_data, thefile)

add_ebps(ACTIVE_LAYER_SUBDIR)
add_ebps(OUTER_LAYER_SUBDIR)