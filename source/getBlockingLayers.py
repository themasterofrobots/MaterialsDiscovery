# -*- coding: utf-8 -*-
import math
import multiprocessing
import json
import os.path
import ast
from pymatgen.core import periodic_table
from pymatgen.matproj.rest import MPRester, MPRestError
from pymatgen.analysis import bond_valence

NUM_PROCESSES = multiprocessing.cpu_count()-1


#Output Settings
PRINT_DEBUG = True
PRINT_INFO = True
SAVE_INFO = True

#These constants currently from https://www.materialsproject.org/wiki/index.php/Calculations_Manual#Band_Structure_and_DOS_computations_details
CORRECTION_SLOPE = 1.6
CORRECTION_OFFSET = 0

EBP_CORRECTION_SLOPE = 1.3403
EBP_CORRECTION_Y_INT = .3968

BAND_GAP_LB = 0.9
BAND_GAP_UB = 7

SINGLE_SUBSCRIPT = 8
TOTAL_SUBSCRIPT = 18

BAND_GAP_INCREMENT = .1
BAND_GAP_UB_TEMP = BAND_GAP_LB + BAND_GAP_INCREMENT

MAX_DEFINED_ELEMENT = 103
MAX_ALLOWED_ELEMENT = 82

FILENAME = 'blockingLayers_binary'

# This exception is an empty exception to handle no material at a given ID
class NoMaterialAtIDError(Exception):
	pass
	#Don't need anything here, just define it as an exception

def calcEBP(bandstruct, name, N_VB, N_CB, PRINT_DEBUG):
	bs = bandstruct
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

	ebp = (vb_en + cb_en)/(2*len(bs.kpoints))
	if not(bs.is_metal()):
		ebp -= (bs.get_vbm())['energy']
	else:
		print 'Warning, material is a metal! No VBM offset applied, no file written'
	if PRINT_DEBUG:
		print 'The branch point energy of ' + name + ' is ' + str(ebp) + ' eV'
	ebp *= EBP_CORRECTION_SLOPE
	ebp += EBP_CORRECTION_Y_INT
	return ebp

def get_materials(filter):
	# Gets materials with specified filter, also calculates branch point enegery, and writes JSON data to list
	with MPRester() as mp:
		query_batch = mp.query(ast.literal_eval(filter), ["pretty_formula", "material_id","bv_structure", "has_bandstructure", "anonymous_formula", "e_above_hull"])
		number_of_compounds = len(query_batch)

		global pool
		results = pool.map(get_material, query_batch)
		allData.extend([elem for elem in results if elem is not None])

def get_material(query_info):
	try:
		#We copy the internals of the get_bandstructure_by_material_id method
		#to handle missing material ids
		with MPRester() as mp:
			if query_info["has_bandstructure"] == False:
				raise NoMaterialAtIDError

			#check if material is stable
			if query_info["e_above_hull"] > 0:
				return

			formula_numbers = query_info["anonymous_formula"]
			first_elem_num = formula_numbers["A"] if "A" in formula_numbers else 0
			second_elem_num = formula_numbers["B"] if "B" in formula_numbers else 0
			third_elem_num = formula_numbers["C"] if "C" in formula_numbers else 0 
			if (first_elem_num > SINGLE_SUBSCRIPT or second_elem_num > SINGLE_SUBSCRIPT 
				or third_elem_num > SINGLE_SUBSCRIPT or first_elem_num + second_elem_num + third_elem_num > TOTAL_SUBSCRIPT):
				return

			name = query_info["pretty_formula"]
			if ('Cd' in name) or ('Hg' in name) or ('Pb' in name) or ('Cr' in name): #check if material has Cd
				if PRINT_DEBUG:
					print "Material with index " + query_info["material_id"] + " contains prohibited material!"
				return

			banddata = mp.get_data(query_info['material_id'], prop="bandstructure")
			bandstruct = banddata[0]["bandstructure"]
				
			bandinfo = bandstruct.get_band_gap()
			bandgap = CORRECTION_SLOPE*bandinfo["energy"]+CORRECTION_OFFSET
			
			if PRINT_INFO:
				print "For " + name + ", we have " + ("a direct" if bandinfo["direct"] else "an indirect") + " band gap of " + str(bandgap) + " eV."

			if (bandstruct.is_metal()):
				return

			# Save information
			oxidizedstates = query_info['bv_structure'].get_primitive_structure()
			total_e = 0
			if(type(oxidizedstates.species[0]) is not periodic_table.Specie):
				return

			#TODO: Investigate how specie_and_occu differs?
			#For now, assume we don't don't
			for specie in oxidizedstates.species:
				if specie.Z > MAX_ALLOWED_ELEMENT:
					return
				num_electrons = specie.Z+specie.oxi_state
				#skip if beyond known elements. This actually works for ions with
				#no electrons (ie H^+ or He^2+)
				if(num_electrons <= 0 or num_electrons > MAX_DEFINED_ELEMENT):
					continue
				#This is a hack of sorts. We get the element that has the same number 
				#of electrons as our ion, and use that as the electronic structure
				orbitals = periodic_table.get_el_sp(int(num_electrons)).full_electronic_structure
				if len(orbitals) <= 1:
					continue
				index = -1
				max_orb_num = orbitals[index][0]
				#This will exclude d and higher electrons because they will never be the top
				#principal quantum number. We could also grab this directly to be safe
				while orbitals[index][0] == max_orb_num:
					total_e += orbitals[index][2]
					index -= 1

			N_VB = total_e/4
			N_CB = total_e/8
			# Calculate branch point energy (ebp) and correct
			ebp = calcEBP(bandstruct, name, N_VB, N_CB, PRINT_DEBUG)
			return {'name':name, 'vbm': -ebp, 'cbm': bandgap-ebp, 'bandinfo': bandinfo}

	except NoMaterialAtIDError:
		if PRINT_DEBUG:
			print "Material with index " + query_info["material_id"] + " does not exist!"

	except MPRestError:
		if PRINT_DEBUG:
			print "No band struct or density of states at index " + query_info["material_id"] + "!"



allData = []
pool = multiprocessing.Pool(NUM_PROCESSES)

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