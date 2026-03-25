NaCl_PP_LR Model
----------------

The file NaCl_PP_LR.mdl can be used with the TENSOAP-FAST code (https://github.com/dilkins/TENSOAP-FAST). Once this code is compiled, the model is applied by using:

/PATH/TO/TENSOAP-FAST/bin/sagpr_apply -f FRAMES.xyz -m /PATH/TO/NaCl_PP_LR.mdl -o prediction_L1.txt

where FRAMES.xyz is an extended xyz file containing the frames for which the polarization is desired. Note that the output will be the lambda=1 spherical component of the polarization *per atom*, so needs to be multiplied by the number of atoms to give the polarization. Alternatively,

/PATH/TO/TENSOAP-FAST/bin/sagpr_apply_process -f FRAMES.xyz -m /PATH/TO/NaCl_PP_LR.mdl -o prediction_L1.txt

gives the lambda=1 spherical component of the total polarization. To convert this to the Cartesian polarization, the following command can be used:

/PATH/TO/TENSOAP-FAST/tools/spherical_to_cartesian.sh 1 prediction_L1.txt prediction_cartesian.txt

this creates a file, prediction_cartesian.txt, which contains the Cartesian polarization (if sagpr_apply_process was used) or polarization per atom (if sagpr_apply was used).
