NaCl_W_SR model
----------------

The model provided here can be applied using the TENSOAP code (https://github.com/dilkins/TENSOAP) using the following steps (where FRAMES.xyz is the file for which predictions are desired):

source /PATH/TO/TENSOAP/env.sh


for lm in 0 1;do
	sagpr_get_PS -f FRAMES.xyz -lm ${lm} -o PS${lm}_out -sf PS${lm} -sg 0.20 -rc 6.0 -rs 1 2.0 5 -c H O Na Cl -s H O Na Cl -p -l 4 -n 6
	get_atomic_power_spectrum.py -p PS${lm}_out.npy -f FRAMES.xyz -o PS${lm}_out_atomic -c O
done
sagpr_get_kernel -z 2 -ps PS1_out_atomic_O.npy PS1_atomic_O.npy -ps0 PS0_out_atomic_O.npy PS0_atomic_O.npy -s NONE NONE -o K1
sagpr_prediction -r 1 -k K1.npy -w weights -o prediction

Then, prediction_cartesian.out contains the predicted Delta_i vectors for each atom. These can be assigned to the atoms using the script provided in ths folder, apply_predicted_wannier_centres.py as:

python3 /PATH/TO/apply_predicted_wannier_centres.py -f FRAMES.xyz -w prediction_cartesian.txt -o wannier.xyz

The output file wannier.xyz is an extended xyz file containing the predicted wannier displacements.
