#!/usr/bin/env python3

import numpy as np
from ase.io import read,write
import argparse,sys,os

# INPUT ARGUMENTS.
parser = argparse.ArgumentParser(description="Wannier centres")
parser.add_argument("-f", "--file", required=True, help="Input frame")
parser.add_argument("-w", "--wannier", required=True, help="File with Wannier centres")
parser.add_argument("-o", "--output", required=True, help="Output file")
parser.add_argument("-e", "--elements", default=["O"], nargs="+", help="List of elements to be assigned centres")
parser.add_argument("-q", "--charges", nargs="+", default=["H","1","O","6"], help="List of charges")
args = parser.parse_args()

frame = read(args.file,":")
w_cen = np.loadtxt(args.wannier)

if (np.shape(w_cen)[1]==6):
  w_cen = w_cen[:,3:]

elements_list = [[] for i in range(len(frame))]
for i in range(len(frame)):
  for j in range(len(frame[i])):
    if frame[i][j].symbol in args.elements:
      elements_list[i].append(j)

#print(frame[0].arrays)

for i in range(len(frame)):
  if "wannier_dist" in frame[i].arrays:
    del frame[i].arrays["wannier_dist"]
  if "Particles+Wannier" in frame[i].info:
    del frame[i].info["Particles+Wannier"]
  if "centers." in frame[i].info:
    del frame[i].info["centers."]
  if "Iteration:1_0" in frame[i].info:
    del frame[i].info["Iteration:1_0"]
  if "wannier_polarization" in frame[i].info:
    del frame[i].info["wannier_polarization"]

if(len(np.concatenate(elements_list)) != len(w_cen)):
  print("ERROR: list of elements to which centres are assigned has a different size to that of data!")
  print(len(np.concatenate(elements_list)),len(w_cen))
  sys.exit(0)

#print(np.shape(w_cen))
k = -1

#print(w_cen[0])

for i in range(len(frame)):
  wannier_dists = np.zeros((len(frame[i]),3),dtype=float)
  for j in range(len(frame[i])):
    if (frame[i][j].symbol in args.elements):
      k += 1
      wannier_dists[j] = w_cen[k]
#      print(j,wannier_dists[j],k,w_cen[k])
#  print(wannier_dists)
#  sys.exit(0)
  frame[i].arrays["wannier_dists"] = wannier_dists

#print(k)

ch_list = {}
for i in range(0,len(frame[0].info["list_of_charges"].split()),2):
  ch_list[frame[0].info["list_of_charges"].split()[i]] = int(frame[0].info["list_of_charges"].split()[i+1])

for i in range(len(frame)):
  qpol = frame[i].get_cell().T * 4.8032047
  pol = np.sum([frame[i][j].position*ch_list[frame[i][j].symbol] for j in range(len(frame[i]))],axis=0) * 4.8032047
  for j in range(len(frame[i])):
    if (frame[i][j].symbol in args.elements):
      pol -= 4.8032047 * 2 * frame[i].info["n_wannier"] * (frame[i][j].position + frame[i].arrays["wannier_dists"][j])
  pol -= np.dot(qpol,np.round(np.dot(np.linalg.inv(qpol),pol),0))
  frame[i].info["predicted_wannier_polarization"] = pol


#check_pol = np.sum([xyz[i].position*ch_list[xyz[i].symbol] for i in range(len(xyz))],axis=0) * 4.8032047

#wannier_vectors = xyz.arrays["wannier_dist"]

#for i in range(len(xyz)):
#  if (xyz[i].symbol in xyz.info["elements_with_centres"].split()):
#    check_pol -= 4.8032047 * 2 * xyz.info["n_wannier"] *(xyz[i].position + wannier_vectors[i])

#check_pol -= np.dot(qpol,np.round(np.dot(np.linalg.inv(qpol),check_pol),0))

write(args.output,frame)
