# -*- coding: utf-8 -*-
from pymatgen.matproj.rest import MPRester, MPRestError
import json
import os.path

#Output Settings
PRINT_DEBUG = True
PRINT_INFO = True
SAVE_INFO = True
SUBFOLDER = 'blockingLayers/'

#These constants currently from https://www.materialsproject.org/wiki/index.php/Calculations_Manual#Band_Structure_and_DOS_computations_details
CORRECTION_SLOPE = 1.6
CORRECTION_OFFSET = 0

#Search Band Gaps (eV)
BAND_GAP_LB = 5.2 #original = 0.8eV 
BAND_GAP_UB = 7

#Test with just Si
# i = 149

#This exception is an empty exception to handle no material at a given ID
class NoMaterialAtIDError(Exception):
	pass
	#Don't need anything here, just define it as an exception

with MPRester() as mp:
	bandgap_filter = mp.query({"band_gap": {"$gt": BAND_GAP_LB, "$lt": BAND_GAP_UB}}, ["pretty_formula", "band_gap", "material_id"])
	# filters out materials with bandgaps outside desired area, pulls formula and materials_id
	number_of_compounds = len(bandgap_filter)
	
	for i in range(0,number_of_compounds-1):
		try:
			#bandstruct = mp.get_bandstructure_by_material_id("mp-" + str(i))
			#We copy the internals of the get_bandstructure_by_material_id method
			#to handle missing material ids
			banddata = mp.get_data(bandgap_filter[i]['material_id'], prop="bandstructure")
			if len(banddata) == 0:
				raise NoMaterialAtIDError
			bandstruct = banddata[0]["bandstructure"]
			#print bandstruct.efermi

			#dos = mp.get_dos_by_material_id("mp-1")

			name = bandgap_filter[i]["pretty_formula"]
			if ('Cd' in name) or ('Hg' in name) or ('Pb' in name) or ('Cr' in name): #check if material has Cd
				if PRINT_DEBUG:
					print "Material with index " + bandgap_filter[i]["material_id"] + " contains prohibited material!"
				continue
			elif sum(1 for numCaps in name if numCaps.isupper()) > 3: #checks number of capital letters in formula
				if PRINT_DEBUG:
					print "Material with index " + bandgap_filter[i]["material_id"] + " is greater than ternary!"
				continue
			formendata = mp.get_data(bandgap_filter[i]['material_id'], prop="formation_energy_per_atom")
			formen = formendata[0]["formation_energy_per_atom"]

			bandinfo = bandstruct.get_band_gap()
			bandgap = CORRECTION_SLOPE*bandinfo["energy"]+CORRECTION_OFFSET
			if PRINT_INFO:
				print "For " + name + ", we have " + ("a direct" if bandinfo["direct"] else "an indirect") + " band gap of " + str(bandgap) + " eV."
			if SAVE_INFO:
					json_filename = SUBFOLDER+name+'_structure.json'
					try:
						with open(json_filename, 'r') as rf:
							prev_json_data = json.load(rf)
							#overwrite if more stable, formation energy is lower
							if PRINT_DEBUG:
								print "Existing file for this material"
							if prev_json_data["formation_energy"] > formen:
								with open(json_filename, 'w') as wf:
									json.dump({'name':name, 'formation_energy':formen, 'bandstructure':bandstruct.as_dict()}, wf)
								if PRINT_DEBUG:
									print "Overwrote existing material file, new material has lower formation energy"
							else:
								if PRINT_DEBUG:
									print "Left existing material file alone, it has lower formation energy"
					except (IOError, ValueError):
						#if file doesn't exist or is malformed, just overwrite it
						if PRINT_DEBUG:
							print "Malformed or no current saved material found, overwritting"
						with open(json_filename, 'w') as wf:
							json.dump({'name':name, 'formation_energy':formen, 'bandstructure':bandstruct.as_dict()}, wf)
		except NoMaterialAtIDError:
			if PRINT_DEBUG:
				print "Material with index " + bandgap_filter[i]["material_id"] + " does not exist!"

		except MPRestError:
			if PRINT_DEBUG:
				print "No band struct or density of states at index " + bandgap_filter[i]["material_id"] + "!"




