models
======

This folder contains the models needed to get the pre-processed and Wannier polarizabilities. To use either of these models, the TENSOAP-FAST code is needed (https://github.com/dilkins/TENSOAP-FAST):

Pre-processed polarizations
---------------------------

To apply the preprocessed polarization model to a file "frames.xyz" in extended-xyz format, use the command:

/PATH/TO/TENSOAP-FAST/bin/sagpr_apply -f frames.xyz -m /path/to/models/unwrap/unwrap.mdl -v -o MU_UNWRAP.out

where unwrap.mdl is in the subfolder "unwrap" here. This will give a file, MU_UNWRAP.out, which has the polarization of each frame PER ATOM, in the form of an L=1 Cartesian tensor (i.e., Py/natom,Pz/natom,Px/natom).

Wannier-centre polarizations
----------------------------

We provide models for the Wannier displacements; these must be post-processed to give the polarization. First, use:

/PATH/TO/TENSOAP-FAST/bin/sagpr_apply -f frames.xyz -m /path/to/models/wannier/wannier_{C,O}.mdl -a -o MU_WANNIER_{C,O}.out

then run /path/to/models/wannier/get_output_with_predicted_centres.sh -- this will give a file called "output_with_predicted_centres.xyz" which is an extended-xyz file that contains the predicted Wannier displacements.

To get the Wannier-centre polarizations run:

python3 /path/to/models/wannier/polarization_from_wannier_centres.py -f output_with_predicted_centres.xyz -q H 1 C 4 O 6 > POL_WANNIER.out

The file "POL_WANNIER.out" will have one line per frame, with the polarization calculated from the positions of Wannier centres.
