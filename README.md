# MaterialsDiscovery
==========

MaterialsDiscovery is a set of scripts written by Scott Biesboer, Brian Cabinian, Andrew Erwin, and Hongyi Michael Wu under the advisement of Professors Andre Schleife and Moonsub Shim at the University of Illinois. MaterialsDiscovery is focused on the data mining of electronic structure information of inorganic compounds from the Materials Project Database in order to identify appropriate heavy-metal-free materials to be used in quantum dot light emitting diode (LED) heterostructures. MaterialsDiscovery is a method for the identification of novel HTL, ETL, and emissive layers and using these materials to design QD visible-light LED heterostructures with potentially high EQEs.

Materials Discovery requires users to have Python and the PyMatGen library (http://pymatgen.org/) to run. Users must also register on https://www.materialsproject.org/ and obtain a API Key, which should be set as the MAPI_KEY environmental variable.

The main files that users should run and edit are CompileMaterialsJSON.py and then find_heterojunctions.py. CompileMaterialsJSON.py will fetch the HTL, ETL, and emissive layers from Materials Project. Then find_heterojunction.py will construct heterojunctions out of these layers, and rank them by injection matching.
The rest of the source files are just for testing and reference, as we built up the script in parts. No guarantee is made to their readability or functionality.

Note that there is a long list of features that should be added, and we encourage people to do so.
1. Use BoltzTraP (Boltzmann Transport Properties) function within Materials Project to calculate effective mass, obtain conductivity
2. Run full band structure calculations on remaining ~2000 compounds to find effective mass, conductivity, mobility, etc.
3. Consider anode and cathode materials
4. Physically manufacture a device to obtain real efficiency data
5. Check that at least one of HTL and ETL is transparent
6. Explore organic materials and relation to branch point energy


