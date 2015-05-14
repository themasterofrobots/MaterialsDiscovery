# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 13:35:07 2015

@author: cabinian
"""
import pymatgen as mg

with mg.MPRester() as m:
    
    #Bandstructure for material id
    bandstructure = m.get_bandstructure_by_material_id("mp-149")
    
#obviously not anything interesting.. don't know how to get stuff out of the bandstructure object
print(bandstructure)