#!/usr/bin/env python3

import numpy as np
from ase.io import read,write
import argparse,sys,os

def MIC(vec,cell,icell):
  # Apply minimum image convention to a vector to account for
  # periodic boundary conditions
  ivc = np.dot(icell,vec)
  rvc = np.round(ivc)
  cvc = np.dot(cell,rvc)
  ovc = vec - cvc
  return ovc

# INPUT ARGUMENTS.
parser = argparse.ArgumentParser(description="Wannier centres")
parser.add_argument("-f", "--file", required=True, help="Input frame")
parser.add_argument("-c", "--cell", required=True, help="File with cell data")
parser.add_argument("-o", "--output", default="output.xyz", help="Output file")
parser.add_argument("-e", "--elements", default=["O"], nargs="+", help="List of elements to be assigned centres")
parser.add_argument("-q", "--charges", nargs="+", default=["H","1","C","4","O","6"], help="List of charges")
args = parser.parse_args()

frame = read(args.file)
cl    = read(args.cell).get_cell().T
ic    = np.linalg.inv(cl)

list_of_elements = []
elements = []
wannier  = []

for i in range(len(frame)):
  if (frame[i].symbol in args.elements):
    list_of_elements.append(i)
    elements.append(frame[i])
  elif (frame[i].symbol == 'X'):
    wannier.append(frame[i])

charges_dict = {'X':-2.0}
for i in range(0,len(args.charges),2):
  charges_dict[args.charges[i]] = float(args.charges[i+1])

# For each Wannier centre, go through and assign it to the closest atom
list_of_centres = [[] for i in range(len(elements))]
for i in range(len(wannier)):
  centre = wannier[i]
  dists = [np.linalg.norm(MIC(centre.position - elements[i].position,cl,ic)) for i in range(len(elements))]
  element = np.argmin(dists)
  list_of_centres[element].append(i)

## Check whether each atom has been assigned the same number of Wannier centres
#for i in range(len(list_of_centres)):
#  if (len(list_of_centres[i]) != len(wannier)/len(elements)):
#    print("ERROR: element",i,"has been assigned a different number of centres than the rest!")
#    print(list_of_centres[i])
#    print(list_of_centres)
#    sys.exit(0)

# Get number of Wannier centres assigned to each atom
n_centres = np.zeros(len(frame),dtype=int)
#print(list_of_elements)
for i in range(len(list_of_elements)):
  n_centres[list_of_elements[i]] = len(list_of_centres[i])
#for i in range(len(frame)):
#  if (frame[i].symbol in args.elements:
#print(len(list_of_centres))

frame.arrays["n_centres"] = n_centres

# Having passed this check, get total vector offset of centres from each atom
diff_vector = [None for i in range(len(list_of_centres))]
for i in range(len(list_of_centres)):
  diff_vector[i] = MIC(np.sum([wannier[list_of_centres[i][j]].position for j in range(len(list_of_centres[i]))],axis=0) - elements[i].position,cl,ic)

# Calculate the Wannier polarization as a check
pol = np.sum([frame[i].position*charges_dict[frame[i].symbol] for i in range(len(frame))],axis=0) * 4.8032047

# Get the quantum of polarization and apply PBCs
qpol = cl * 4.8032047
pol -= np.dot(qpol,np.round(np.dot(np.linalg.inv(qpol),pol),0))

# Go through and put this information into the array
all_wannier = np.zeros((len(frame),3),dtype=float)
for i in range(len(frame)):
  if (i in list_of_elements):
    all_wannier[i] = diff_vector[list_of_elements.index(i)]
frame.arrays["wannier_dist"] = all_wannier
frame.set_cell(cl.T)

# Get rid of Wannier centres
for i in range(len(frame)-1,0,-1):
  if (frame[i].symbol=='X'):
    frame.pop(i)

# Print file
#frame.info["n_wannier"] = len(wannier)/len(elements)
frame.info["elements_with_centres"] = ' '.join(args.elements)
frame.info["list_of_charges"] = ' '.join(args.charges)
frame.info["wannier_polarization"] = pol
frame.set_pbc(True)
write(args.output,frame)

# As a final check, read in the output file and make sure we can recalculate the polarization from it
xyz = read(args.output)
ch_list = {}
for i in range(0,len(xyz.info["list_of_charges"].split()),2):
  ch_list[xyz.info["list_of_charges"].split()[i]] = int(xyz.info["list_of_charges"].split()[i+1])

check_pol = np.sum([xyz[i].position*ch_list[xyz[i].symbol] for i in range(len(xyz))],axis=0) * 4.8032047

wannier_vectors = xyz.arrays["wannier_dist"]

for i in range(len(xyz)):
  if (xyz[i].symbol in xyz.info["elements_with_centres"].split()):
    check_pol -= 4.8032047 * 2 *(xyz[i].position + wannier_vectors[i])

check_pol -= np.dot(qpol,np.round(np.dot(np.linalg.inv(qpol),check_pol),0))

if (np.linalg.norm(check_pol - xyz.info["wannier_polarization"]) > 1e-5):
  print("ERROR: mismatch in polarizations:")
  print(check_pol,xyz.info["wannier_polarization"])
  print(np.linalg.norm(check_pol-xyz.info["wannier_polarization"]))
  sys.exit(0)
