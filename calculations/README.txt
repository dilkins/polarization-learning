cp2k calculations
-----------------

The file cp2k.in is used to carry out cp2k calculations, outputting the positions of Wannier centres as well as the total polarization computed using these centres. The six variables at the beginning of the input file must be appropriately set to give the correct simulation box. Calculations were carried out with CP2K version 4.1.

Each of the three subfolders, 01_bulk_water, 02_interfacial_water and 03_concentrated_nacl, contains the frames on which calculations were run. For each frame, the polarization calculated from Wannier centres ("wannier_polarization") is given, along with the Wannier deltas for each oxygen atom (columns 5-7 for each atom entry).
