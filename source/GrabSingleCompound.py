# -*- coding: utf-8 -*-
from pymatgen.matproj.rest import MPRester, MPRestError
import pymatgen.analysis.bond_valence
from pymatgen.core import periodic_table
import json
import os.path

#Output Settings
PRINT_DEBUG = True
PRINT_INFO = True
SAVE_INFO = True

#These constants currently from https://www.materialsproject.org/wiki/index.php/Calculations_Manual#Band_Structure_and_DOS_computations_details
CORRECTION_SLOPE = 1.6
CORRECTION_OFFSET = 0

#ZnO
#NAME='mp-2133'

#AlN
#NAME='mp-661'

#GaN-wz
#NAME='mp-804'

#GaN-zb
#NAME='mp-830'

#MgO
#NAME='mp-1265'

#CdO
#NAME='mp-1132'

#In2O3
#NAME='mp-22323'

NAME='mp-647209'

MAX_DEFINED_ELEMENT = 103

#This exception is an empty exception to handle no material at a given ID
class NoMaterialAtIDError(Exception):
	pass
	#Don't need anything here, just define it as an exception

with MPRester() as mp:
	try:
		#bandstruct = mp.get_bandstructure_by_material_id("mp-" + str(i))
		#We copy the internals of the get_bandstructure_by_material_id method
		#to handle missing material ids
		banddata = mp.get_data(NAME, prop="bandstructure")
		if len(banddata) == 0:
			raise NoMaterialAtIDError

		oxidizedstates = mp.query(NAME, ["bv_structure"])[0]['bv_structure'].get_primitive_structure()

		total_e = 0

		#TODO: Investigate how specie_and_occu differs?
		#For now, assume we don't don't
		for specie in oxidizedstates.species:
			num_electrons = specie.Z+specie.oxi_state
			#skip if beyond known elements. This actually works for ions with
			#no electrons (ie H^+ or He^2+)
			if(num_electrons <= 0 or num_electrons > MAX_DEFINED_ELEMENT):
				continue

			#This is a hack of sorts. We get the element that has the same number 
			#of electrons as our ion, and use that as the electronic structure
			orbitals = periodic_table.get_el_sp(num_electrons).full_electronic_structure
			index = -1
			max_orb_num = orbitals[index][0]
			#This will exclude d and higher electrons because they will never be the top
			#principal quantum number. We could also grab this directly to be safe
			while orbitals[index][0] == max_orb_num:
				total_e += orbitals[index][2]
				index -= 1

		N_VB = total_e/4
		N_CB = total_e/8
		if(PRINT_DEBUG):
			print "Total e = " + str(total_e) + ", N_VB = " + str(N_VB) + ", N_CB = " + str(N_CB)


		formendata = mp.get_data(NAME, prop="formation_energy_per_atom")
		for i in range(len(banddata)):
			bandstruct = banddata[i]["bandstructure"]
			#print bandstruct.efermi

			#dos = mp.get_dos_by_material_id("mp-1")

			formen = formendata[i]["formation_energy_per_atom"]

			bandinfo = bandstruct.get_band_gap()
			bandgap = CORRECTION_SLOPE*bandinfo["energy"]+CORRECTION_OFFSET
			if PRINT_INFO:
				print "For " + NAME + ", we have " + ("a direct" if bandinfo["direct"] else "an indirect") + " band gap of " + str(bandgap) + " eV."
			if SAVE_INFO:
				json_filename = NAME+'_structure.json'
				try:
					with open(json_filename, 'r') as rf:
						prev_json_data = json.load(rf)
						#overwrite if more stable, formation energy is lower
						if PRINT_DEBUG:
							print "Existing file for this material"
						if prev_json_data["formation_energy"] > formen:
							with open(json_filename, 'w') as wf:
								json.dump({'name':NAME, 'formation_energy':formen, 'bandstructure':bandstruct.as_dict()}, wf)
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
						json.dump({'name':NAME, 'formation_energy':formen, 'bandstructure':bandstruct.as_dict()}, wf)

	except NoMaterialAtIDError:
		if PRINT_DEBUG:
			print "Material " + NAME + " does not exist!"

	except MPRestError:
		if PRINT_DEBUG:
			print "No band struct or density of states for " + NAME + "!"




