inputs
======

This folder contains one subfolder, with the inputs needed to run CP2K and scripts to post-process it. CP2K must be run in a folder containing an xyz file, "frame.xyz", and the input files in the "cp2k" subfolder. The output should go into a file called wannier.out.

Running the script "process.sh" in the folder "cp2k-postproc" will produce a file called "output.xyz", which is in extended-xyz format and will contain both the raw polarizability and the Wannier displacements (described in the main text).
