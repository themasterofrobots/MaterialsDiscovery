from __future__ import division
import math
import json
import os.path
from pymatgen.electronic_structure.bandstructure import BandStructure

#These constants currently from https://www.materialsproject.org/wiki/index.php/Calculations_Manual#Band_Structure_and_DOS_computations_details
CORRECTION_SLOPE = 1.6
CORRECTION_OFFSET = 0

EBP_CORRECTION = 1.5

PRINT_DEBUG = False

ACTIVE_LAYER_SUBDIR = 'activeLayer'
BLOCKING_LAYERS_SUBDIR = 'blockingLayers'

activelayers = []
blockinglayers = []

def add_materials(dirname, data_list):
	for dir_entry in os.listdir(dirname):
		add_material(dirname + '/' + dir_entry, data_list)

def add_material(dir_entry, data_list):
	with open(dir_entry, 'r') as thefile:
			json_data = json.load(thefile)
			if 'ebp' not in json_data:
				print 'Skipped ' + json_data['name'] + ', no ebp'
				return
			bs = BandStructure.from_dict(json_data['bandstructure'])
			bandinfo = bs.get_band_gap()
			bandgap = CORRECTION_SLOPE*bandinfo["energy"]+CORRECTION_OFFSET
			ebp = EBP_CORRECTION*json_data['ebp']
			if PRINT_DEBUG:
				print json_data['name'] + ' has vbm ' + str(-ebp) + " eV and cbm " + str(bandgap-ebp) + ' eV'
			data_list.append({'name': json_data['name'], 'vbm': -ebp, 'cbm': bandgap-ebp})

print 'Started adding active layers'
add_materials(ACTIVE_LAYER_SUBDIR, activelayers)
print 'Finished adding active layers'

with open('activelayers_list.json', 'w') as thefile:
	json.dump(activelayers, thefile)

print 'Started adding blocking layers'
add_materials(BLOCKING_LAYERS_SUBDIR, blockinglayers)
print 'Finished adding blocking layers'

with open('blockinglayers_list.json', 'w') as thefile:
	json.dump(blockinglayers, thefile)