# -*- coding: utf-8 -*-
from __future__ import division
import os
import json
import math
import numpy
from pymatgen.electronic_structure.boltztrap import BoltztrapRunner, BoltztrapAnalyzer
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.electronic_structure.core import Spin
from pymatgen.matproj.rest import MPRester
from pymatgen.core.structure import Structure

#Output Settings
PRINT_DEBUG = True
PRINT_INFO = True
SAVE_INFO = True

"""
1. Read material id from json file
2. Query to get materials project band structure object
	1. Get band structure symm line object
	2. Get structure object from MatProj
	3. Combine the two using from_dict method to get a band structure object that has the 
	proper attributes
3. Run BoltztrapRunner on the band structure object
4. Extract conductivity data using from_files
"""

with MPRester() as mp:
	#banddata = mp.query(criteria={"task_id": "mp-804"}, properties=["band_structure_uniform"])
	banddata = mp.get_data("mp-804", prop="bandstructure")
	bs = banddata[0]["bandstructure"]
	structdata = mp.query("mp-804", ["structure"])
	bs._structure = structdata[0]["structure"]
	#bs = mp.get_bandstructure_by_material_id('mp-804')

	i = 0
	while i < len(bs._kpoints):
		for band in bs._bands[Spin.up]:
			band.pop(i)
		bs._kpoints.pop(i)
		i = i + 1

	existkpoints = []
	i = 0
	while i < len(bs._kpoints):
		if any(numpy.array_equal(bs._kpoints[i].frac_coords, kpt) for kpt in existkpoints):
			for band in bs._bands[Spin.up]:
				band.pop(i)
			bs._kpoints.pop(i)
		else:
			existkpoints.append(bs._kpoints[i].frac_coords)
			i = i + 1


	nelec = 38 

	x = BoltztrapRunner(bs, nelec, dos_type="HISTO", energy_grid=0.005, lpfac=10, type="BOLTZ", band_nb=None, tauref=0, tauexp=0, tauen=0, soc=False)
	#goes to nonetype... band structure objects have a ._structure attribute but not .structure
	outfile = x.run(prev_sigma=None, path_dir=os.getcwd(), convergence=True)

	y = BoltztrapAnalyzer.from_files(outfile)
		#adjusted doping levels in boltztrap.outputtrans file in CoSb3 folder
	EffMass = y.get_eig_average_eff_mass_tensor(temperature=300, doping=1e+18)
