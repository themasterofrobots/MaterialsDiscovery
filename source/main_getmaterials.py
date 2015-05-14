# -*- coding: utf-8 -*-
import math
import multiprocessing
import json
import os.path
import ast
from pymatgen.core import periodic_table
from pymatgen.matproj.rest import MPRester, MPRestError
from pymatgen.analysis import bond_valence
from CompileMaterialsJSON import NoMaterialAtIDError, calcEBP, get_materials, get_material

# Output Settings
PRINT_DEBUG = True
PRINT_INFO = True
SAVE_INFO = True

# ------- Global settings when finding materials for active and blocking layers -------
# These constants currently from https://www.materialsproject.org/wiki/index.php/Calculations_Manual#Band_Structure_and_DOS_computations_details
CORRECTION_SLOPE = 1.6
CORRECTION_OFFSET = 0

# These are the correction factors for the EBP to accomodate for the fact that the
# MaterialProject only has high symmetry paths, but we need the full Brillouin Zone
EBP_CORRECTION_SLOPE = 1.3403
EBP_CORRECTION_Y_INT = .3968

# Force the subscripts to be small to only get simple compounds
global SINGLE_SUBSCRIPT
global TOTAL_SUBSCRIPT
SINGLE_SUBSCRIPT = 8
TOTAL_SUBSCRIPT = 18

#Limit the elements to not radioactive compounds and ions that aren't too charged
global MAX_DEFINED_ELEMENT
global MAX_ALLOWED_ELEMENT
MAX_DEFINED_ELEMENT = 103
MAX_ALLOWED_ELEMENT = 82

# ------- Find active layer materials --------
BAND_GAP_LB = .9
BAND_GAP_UB = 3
CHECKDIRECT = True
FILENAME = 'activeLayers'
activelayerfilters = ['{"elements":{"$in":["N", "P", "As", "Sb", "Bi"], "$all":["Ga"]}, "nelements":{"$lte":2}, "band_gap": {"$gt": '+ str(BAND_GAP_LB) + ', "$lt": ' + str(BAND_GAP_UB) +'}}',
 			   '{"elements":{"$in":["N", "P", "As", "Sb", "Bi"], "$all":["In"]}, "nelements":{"$lte":2}, "band_gap": {"$gt": '+ str(BAND_GAP_LB) + ', "$lt": ' + str(BAND_GAP_UB) +'}}',
 				'{"elements":{"$all":["In", "Ga"]}, "nelements":{"$lte":2}, "band_gap": {"$gt": '+ str(BAND_GAP_LB) + ', "$lt": ' + str(BAND_GAP_UB) +'}}',
 				'{"elements":{"$all":["Cu", "In"]}, "nelements":{"$lte":2}, "band_gap": {"$gt": '+ str(BAND_GAP_LB) + ', "$lt": ' + str(BAND_GAP_UB) +'}}',
 				'{"elements":{"$all":["Ag", "Zn"]}, "nelements":{"$lte":2}, "band_gap": {"$gt": '+ str(BAND_GAP_LB) + ', "$lt": ' + str(BAND_GAP_UB) +'}}',
 				'{"elements":{"$all":["Ag", "Ga"]}, "nelements":{"$lte":2}, "band_gap": {"$gt": '+ str(BAND_GAP_LB) + ', "$lt": ' + str(BAND_GAP_UB) +'}}',
 				'{"elements":{"$in":["O", "S", "Se", "Te"], "$all":["Zn"]}, "nelements":{"$lte":2}, "band_gap": {"$gt": '+ str(BAND_GAP_LB) + ', "$lt": ' + str(BAND_GAP_UB) +'}}']

allData = []

for i_filter in range(len(activelayerfilters)):
	print 'Started using filter ' + str(i_filter+1) + '.'
	get_materials(activelayerfilters[i_filter])

json_filename = FILENAME + '.json'
with open(json_filename, 'w') as wf:
	json.dump(allData, wf)

# ------- Find blocking layer materials --------
BAND_GAP_LB = 0.9
BAND_GAP_UB = 7
BAND_GAP_INCREMENT = .1
BAND_GAP_UB_TEMP = BAND_GAP_LB + BAND_GAP_INCREMENT
CHECKDIRECT = False
FILENAME = 'blockingLayers'

allData = []

for i_filter in range(int((BAND_GAP_UB-BAND_GAP_LB)/BAND_GAP_INCREMENT)):
	print 'Started using filter ' + str(i_filter+1) + '.'
	get_materials('{"elements":{"$nin":["Cd", "Hg", "Cr", "Pb"]}, "nelements":{"$lte":2}, "band_gap": {"$gte": '+ str(BAND_GAP_LB) + ', "$lt": ' + str(BAND_GAP_UB_TEMP) +'}}')
	BAND_GAP_LB += BAND_GAP_INCREMENT
	BAND_GAP_UB_TEMP += BAND_GAP_INCREMENT
	print BAND_GAP_LB
	print BAND_GAP_UB_TEMP

json_filename = FILENAME + '.json'
with open(json_filename, 'w') as wf:
	json.dump(allData, wf)