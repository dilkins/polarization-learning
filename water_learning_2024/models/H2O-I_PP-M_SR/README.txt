H2O-I_PP-M_SR model
-------------------

The model provided here can be applied using the TENSOAP code (https://github.com/dilkins/TENSOAP) using the following steps (where FRAMES.xyz is the file for which predictions are desired):

source /PATH/TO/TENSOAP/env.sh
for lm in 0 1;do
	sagpr_get_PS -f FRAMES.xyz -lm ${lm} -p -sf PS${lm}_train -sg 0.20 -rc 5.0 -c H O -s H O -rs 1 2.0 5-o PS${lm}_out
done
sagpr_get_kernel -z 2 -ps PS1_out.npy PS1_active.npy -ps0 PS0_out.npy PS0_active.npy -s PS1_out_natoms.npy NONE -o K1
sagpr_prediction -r 1 -k K1.npy -w weights -o prediction

The file prediction_cartesian.out contains the predicted polarizations for each frame.
