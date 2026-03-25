#!/usr/bin/env python3

import numpy as np
from ase.io import read,write
import argparse,sys

# INPUT ARGUMENTS.
parser = argparse.ArgumentParser(description="Get polarization from Wannier centres")
parser.add_argument("-f", "--file", required=True, help="Input file")
parser.add_argument("-q", "--charges", nargs="+", required=True, help="List of charges")
args = parser.parse_args()

frames = read(args.file,":")
nfr = len(frames)

charges = {'X':-2.0}
for i in range(0,len(args.charges),2):
  charges[args.charges[i]] = float(args.charges[i+1])

for i in range(nfr):
  frame = frames[i]
  qpol = frame.get_cell().T * 4.8032047
  wannier_dist = frame.arrays["wannier_dist"]
  n_wannier    = frame.arrays["n_wannier"]

  # Get total polarization
  pol = np.sum([frame[j].position*charges[frame[j].symbol] for j in range(len(frame))],axis=0) * 4.8032047
  for j in range(len(frame)):
    pol -= 4.8032047 * 2 * n_wannier[j] * (frame[j].position + wannier_dist[j])

  # Apply PBCs
  pol -= np.dot(qpol,np.round(np.dot(np.linalg.inv(qpol),pol),0))
  print(' '.join(pol.astype(str)))

sys.exit(0)


charges = {'H':1.0,'O':6.0}
n_w     = {'H':0,'O':4}
qq = 2.5417465*1.8897261

for i in range(nfr):
  frame = frames[i]
  qpol  = frame.get_cell().T * qq
  wannier_dist = frame.arrays["wannier_dists"]
  n_wannier    = [n_w[frame[j].symbol] for j in range(len(frame))]



  # Get total polarization
  check_pol = np.sum([frame[j].position*charges[frame[j].symbol] for j in range(len(frame))],axis=0) * 4.8032047
  for j in range(len(frame)):
    check_pol -= 4.8032047 * 2 * n_wannier[j] * (frame[j].position + wannier_dist[j])
  check_pol -= np.dot(qpol,np.round(np.dot(np.linalg.inv(qpol),check_pol),0))

  # Build this up from partial polarizations
  nuclear_pol = np.array( [frame[j].position*charges[frame[j].symbol]*qq for j in range(len(frame))])
  electro_pol = np.array( [-2*qq*n_wannier[j] * (frame[j].position + wannier_dist[j]) for j in range(len(frame))])
  check_pol = np.sum(nuclear_pol,axis=0) + np.sum(electro_pol,axis=0)
  check_pol -= np.dot(qpol,np.round(np.dot(np.linalg.inv(qpol),check_pol),0))

  # Resum to get molecular dipole moments
  total_pol = (nuclear_pol + electro_pol).reshape(int(len(frame)/3),3,3)
  total_pol = np.sum(total_pol,axis=1)
  for j in range(len(total_pol)):
    print(' '.join(total_pol[j].astype(str)))
