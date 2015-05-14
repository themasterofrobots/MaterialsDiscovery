from __future__ import division
from operator import itemgetter
import math
import json
import os.path
from pymatgen.electronic_structure.bandstructure import BandStructure

import multiprocessing

NUM_PROCESSES = multiprocessing.cpu_count()-1

#These constants currently from https://www.materialsproject.org/wiki/index.php/Calculations_Manual#Band_Structure_and_DOS_computations_details
CORRECTION_SLOPE = 1.6
CORRECTION_OFFSET = 0

EBP_CORRECTION = 1.5

HOLE_BLOCKING_LTHRESHOLD = 2.5
ELECTRON_BLOCKING_LTHRESHOLD = 2.5
HOLE_INJECTION_LTHRESHOLD = -0.5
ELECTRON_INJECTION_LTHRESHOLD = -0.5
HOLE_BLOCKING_UTHRESHOLD = 10
ELECTRON_BLOCKING_UTHRESHOLD = 10
HOLE_INJECTION_UTHRESHOLD = 1
ELECTRON_INJECTION_UTHRESHOLD = 1
MATCHING_THRESHOLD = 0.075
PRINT_DEBUG = True
PRINT_UNIQUE = True
SAVE_UNIQUE = True

NUM_COMBOS_PRINT = 50
activeLayersFile = 'activelayers.json'
blockingLayersFile = 'blockinglayers.json'
candidatesFile = 'combos.txt'


def findhetero(al):
	global candidate_count
	global blockinglayers
	listofAL = []
	listofHTL = []
	listofETL = []
	listofcombos = []
	if PRINT_DEBUG:
		print 'Checking active layer ' + al['name']
	for etl in blockinglayers:
		if (al['vbm'] - etl['vbm'] > HOLE_BLOCKING_LTHRESHOLD and 
			al['vbm'] - etl['vbm'] < HOLE_BLOCKING_UTHRESHOLD and
			etl['cbm'] - al['cbm'] > ELECTRON_INJECTION_LTHRESHOLD and
			etl['cbm'] - al['cbm'] < ELECTRON_INJECTION_UTHRESHOLD ):
			for htl in blockinglayers:
				injection_matching = abs(abs(etl['cbm'] - al['cbm']) - abs(al['vbm'] - htl['vbm']))
				if (htl['cbm'] - al['cbm'] > ELECTRON_BLOCKING_LTHRESHOLD and
					htl['cbm'] - al['cbm'] < ELECTRON_BLOCKING_UTHRESHOLD and 
					al['vbm'] - htl['vbm'] > HOLE_INJECTION_LTHRESHOLD and
					al['vbm'] - htl['vbm'] < HOLE_INJECTION_UTHRESHOLD and 
					injection_matching < MATCHING_THRESHOLD):
					with candidate_count.get_lock():
						candidate_count.value += 1
						if PRINT_DEBUG:
							print 'Candidate: ' + etl['name'] + ' - ' + al['name'] + ' - ' + htl['name']
							print '      ' + str(etl['cbm']) + ' - ' + str(al['cbm']) + ' - ' + str(htl['cbm'])
							print '      ' + str(etl['vbm']) + ' - ' + str(al['vbm']) + ' - ' + str(htl['vbm'])
						listofAL.append(str(al['name']))
						listofHTL.append(str(htl['name']))
						listofETL.append(str(etl['name']))
						listofcombos.append([injection_matching, str(etl['name']), str(al['name']), str(htl['name'])])
	return [listofAL, listofHTL, listofETL, listofcombos]

activelayers = []
blockinglayers = []


with open(activeLayersFile, 'r') as thefile:
	activelayers = json.load(thefile)

with open(blockingLayersFile, 'r') as thefile:
	blockinglayers = json.load(thefile)

if PRINT_DEBUG:
	print 'We are using ' + str(NUM_PROCESSES) + ' processes.'
	print 'There are ' + str(len(activelayers)) + ' active layers.'
	print 'There are ' + str(len(blockinglayers)) + ' blocking layers.'

print 'Started finding combinations'

candidate_count = multiprocessing.Value('i', 0)
p = multiprocessing.Pool(NUM_PROCESSES)
matslist = p.map(findhetero, activelayers)

#for al in activelayers:
#	findhetero(al)

listofAL = []
listofHTL = []
listofETL = []
listofcombos = []

for sublist in matslist:
	listofAL.extend(sublist[0])
	listofHTL.extend(sublist[1])
	listofETL.extend(sublist[2])
	listofcombos.extend(sublist[3])


#remove duplicates
listofAL = list(set(listofAL))
listofHTL = list(set(listofHTL))
listofETL = list(set(listofETL))

if PRINT_UNIQUE:
	print '\n====================\n'
 	print 'Unique active layers:'
 	print listofAL
 	print 'Unique hole transport layers:'
 	print listofHTL
 	print 'Unique electron transport layers:'
 	print listofETL

print '\n====================\n'
print 'Final Count: ' + str(candidate_count.value) + ' valid heterostructures'

sortedcombos = sorted(listofcombos, key=itemgetter(0))

with open(candidatesFile, 'w') as wf:
	wf.write('Final Count: ' + str(candidate_count.value) + ' valid heterostructures\n')
	wf.write('Top candidates:\n')
	i = 0
	while i < NUM_COMBOS_PRINT and i < len(sortedcombos):
		print ('Candidate: ' + sortedcombos[i][1] + ' - ' + sortedcombos[i][2] + ' - ' +
				sortedcombos[i][3] + ' has injection mismatch of ' + str(sortedcombos[i][0]))
		wf.write('Candidate: ' + sortedcombos[i][1] + ' - ' + sortedcombos[i][2] + ' - ' +
				sortedcombos[i][3] + ' has injection mismatch of ' + str(sortedcombos[i][0]) + '\n')
		i = i + 1

if SAVE_UNIQUE:
	with open('unique_als.json', 'w') as wf:
		json.dump(listofAL, wf)
	with open('unique_htls.json', 'w') as wf:
		json.dump(listofHTL, wf)
	with open('unique_etls.json', 'w') as wf:
		json.dump(listofETL, wf)
