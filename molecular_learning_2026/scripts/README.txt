scripts
=======

This folder contains a script for unwrapping the polarization in a data-driven way, as described in the manuscript. It requires the TENSOAP code () to run, and before starting this script /path/to/TENSOAP/env.sh should have been run to put all key TENSOAP scripts into the PATH variable.

As written, it takes in a file, "train_800.xyz", which contains 800 frames in extended-xyz format, with an L=1 SA-GPR kernel K1_NM_800.npy (between the 800 training frames and the members of an active set). The output, train_UNWRAPPED.xyz, will have the polarizations (which should be denoted "mu" in the header lines for each frame) unwrapped to lie on the same branch.
