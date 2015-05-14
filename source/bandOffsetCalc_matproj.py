from __future__ import division
import math
import json
import os.path
from pymatgen.electronic_structure.bandstructure import BandStructure

#N_CB = 2
#N_VB = 4
#filename = 'ZnO_structure.json'

#N_CB = 2
#N_VB = 4
#filename = 'wz-GaN_structure.json'

#N_CB = 1
#N_VB = 2
#filename = 'MgO_structure.json'

#N_CB = 2
#N_VB = 4
#filename = 'AlN_structure.json'

#N_CB = 1
#N_VB = 2
#filename = 'zb-GaN_structure.json'

N_CB = 1
N_VB = 2
filename = 'CdO_structure.json'

with open(filename, 'r') as thefile:
	json_data = json.load(thefile)
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
	else:
		print "Warning, material is a metal! No VBM offset applied."
	print total_en
