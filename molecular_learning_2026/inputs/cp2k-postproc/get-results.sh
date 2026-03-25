#!/bin/bash

grep -A1 "Dipole moment \[Debye\]" wannier.out | awk '/X=/{printf "%.6f %.6f %.6f\n",$2,$4,$6}' > mu.out
python3 -c $'from ase.io import read,write;import numpy as np;xyz = read("wannier.xyz");fr=read("frame.xyz");xyz.info["description"]=fr.info["description"];xyz.info["frame"]=fr.info["frame"];xyz.info["mu"]=np.loadtxt("mu.out");write("output.xyz",xyz)'
cat output.xyz | sed "s/Particles+Wannier=T centers.=T Iteration:1_0=T//" > tp
mv tp output.xyz
