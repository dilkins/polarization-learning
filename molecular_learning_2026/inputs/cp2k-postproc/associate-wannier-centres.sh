#!/bin/bash

RUNDIR=$(dirname $(readlink -f "$0"))

python3 ${RUNDIR}/associate_wannier_centres_unrestricted.py -f wannier-HOMO_centers_s1-1_0.xyz -c frame.xyz -o wannier.xyz -e C O
