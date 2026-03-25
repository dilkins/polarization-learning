#!/usr/bin/env python3

import numpy as np
from ase.io import read,write
import diptest,subprocess,random,os,warnings
import scipy.linalg
import unwrap_utils

warnings.filterwarnings("ignore")

# Keep going until the distribution of residuals is multimodal
dip = 0.0

k = np.load("K1_NM_800.npy")
nfr = len(k)

print("GETTING POINTS ON THE MAIN BRANCH")
j = 0
while (dip<0.05):

  # Randomize frames
  xyz = read("train_800.xyz",":")

  for i in range(len(xyz)):
    xyz[i].info["frame_num"] = i+1

  write("train_IN.xyz",xyz)

  # Shuffle kernel and data

  idx = [i for i in range(len(xyz))]
  random.shuffle(idx)
  k_out = np.array([k[idx[i]] for i in range(len(idx))])
  xyz_rdm = [xyz[idx[i]] for i in range(len(idx))]
  np.save("K1_NM_RDM.npy",k_out)
  write("train_RDM.xyz",xyz_rdm)

  # Do training with TENSOAP

  subprocess.run(["sagpr_train","-sel","0","1","-r","1","-p","mu","-pr","-perat","-reg","1e-10","-w","weights_1_1e-10","-f","train_RDM.xyz","-sf","K1_NM_RDM.npy","K1_MM.npy"],stdout = subprocess.DEVNULL,stderr = subprocess.DEVNULL)

  # Read in predictions and get residuals
  data = np.loadtxt("prediction_cartesian.txt")[:,:6]
  res = data[:,3:] - data[:,:3]
  # Quantum of polarization
  qpol = [4.8032047 * xyz_rdm[i].get_cell().T for i in range(len(xyz_rdm))]
  # Normalize residuals by QPOL.
  res_norm = [np.dot(np.linalg.inv(qpol[i+1]),res[i]) for i in range(len(res))]
  res = np.array([np.linalg.norm(res_norm[i])**2 for i in range(len(res_norm))])

  dip,pval = diptest.diptest(res)
  print(dip)
  j += 1
  fl = open("scatterplot_" + str(j) + ".txt","w")
  print("#",dip,file=fl)
  for i in range(len(data)):
    print(data[i,0],data[i,3],file=fl)
    print(data[i,1],data[i,4],file=fl)
    print(data[i,2],data[i,5],file=fl)
  fl.close()

# Do training on the main branch
print("TRAINING MODEL ON THE MAIN BRANCH")

frn = np.array([xyz_rdm[i].info["frame_num"] for i in range(len(xyz_rdm))])

# Get list of frames to keep on the main branch
rmx = [0] + [np.max(np.abs(res_norm[i])) for i in range(len(res_norm))]
keep_list = [frn[i]-1 for i in range(len(frn)) if rmx[i]<0.5]
np.savetxt("keep_list.txt",keep_list)

# Get kernel and frames for the main branch
xyz_in = read("train_IN.xyz",":")
write("train_OUT.xyz",[xyz_in[keep_list[i]] for i in range(len(keep_list))])
np.save("K1_NM_OUT.npy",np.array([k[keep_list[i]] for i in range(len(keep_list))]))

# Train model only with points on the main branch
ntr_main = int(len(np.loadtxt("keep_list.txt"))*0.8)
proc1 = subprocess.Popen(["sagpr_train","-rdm",str(ntr_main),"-r","1","-p","mu","-pr","-perat","-reg","1e-5","-w","weights_OUT","-f","train_OUT.xyz","-sf","K1_NM_OUT.npy","K1_MM.npy"],stderr = subprocess.DEVNULL,stdout = subprocess.PIPE)
proc2 = subprocess.Popen(["grep","%"],stdin = proc1.stdout,stdout = subprocess.PIPE,stderr = subprocess.PIPE)
proc1.stdout.close()
out,err = proc2.communicate()
print(out.decode("utf-8").rstrip())

# Get delta-n values
unwrap_list = [frn[i]-1 for i in range(len(frn)) if rmx[i]>0.5]
np.save("K1_NM_UNWRAP.npy",np.array([k[unwrap_list[i]] for i in range(len(unwrap_list))]))

subprocess.run(["sagpr_prediction","-w","weights_OUT","-r","1","-k","K1_NM_UNWRAP.npy","-o","prediction_to_unwrap_peratom"],stdout = subprocess.DEVNULL,stderr = subprocess.DEVNULL)

cartesian = np.loadtxt("prediction_to_unwrap_peratom_cartesian.txt")
output = [cartesian[i]*len(xyz[unwrap_list[i]]) for i in range(len(unwrap_list))]

res = output-np.array([xyz[unwrap_list[i]].info["mu"] for i in range(len(unwrap_list))])

qpol = [4.8032047 * xyz[i].get_cell().T for i in range(len(xyz))]
delta_n = np.round([np.dot(np.linalg.inv(qpol[unwrap_list[i]]),res[i]) for i in range(len(unwrap_list))],0).astype(int)
np.savetxt("delta_n.txt",delta_n)

print("TRAINING FINAL MODEL")

# Apply delta n values to unwrap the data

delta_n = np.loadtxt("delta_n.txt").astype(int)
for i in range(len(xyz)):
  xyz[i].info["mu_unwrapped"] = xyz[i].info["mu"]

for i in range(len(unwrap_list)):
  xyz[unwrap_list[i]].info["mu_unwrapped"] += np.dot(4.8032047*xyz[unwrap_list[i]].get_cell().T,delta_n[i])
write("train_UNWRAPPED.xyz",xyz)

# Train final model

proc1 = subprocess.Popen(["sagpr_train","-rdm","600","-r","1","-p","mu_unwrapped","-pr","-perat","-reg","1e-5","-w","weights_UNWRAPPED","-f","train_UNWRAPPED.xyz","-sf","K1_NM_800.npy","K1_MM.npy"],stderr = subprocess.DEVNULL,stdout = subprocess.PIPE)
proc2 = subprocess.Popen(["grep","%"],stdin = proc1.stdout,stdout = subprocess.PIPE,stderr = subprocess.PIPE)
proc1.stdout.close()
out,err = proc2.communicate()
print(out.decode("utf-8").rstrip())

prediction,mu_test = unwrap_utils.train(k,np.load("K1_MM.npy"),int(0.8*nfr),1e-5,np.array([xyz[i].info["mu_unwrapped"] / len(xyz[i]) for i in range(len(xyz))]))

print(100 * np.sqrt(np.mean((prediction - mu_test)**2)) / np.sqrt(np.mean(mu_test**2)),"%")

k_nm = np.load("K1_NM_800.npy")
k_mm = np.load("K1_MM.npy")
k_nm_train = k_nm[:600]
k_nm_test = k_nm[600:]
nf = len(k_nm)
ne = len(k_nm[0])
nt = len(k_nm_train)
k_nm_train = k_nm_train.transpose(0,2,1,3).reshape((3*nt,3*ne))
k_nm_test  = k_nm_test.transpose(0,2,1,3).reshape((3*(nf-nt),3*ne))
k_mm = k_mm.transpose(0,2,1,3).reshape((3*ne,3*ne))
kmnnm = np.dot(k_nm_train.T,k_nm_train)
kmnnm_reg = kmnnm + 1e-5 * k_mm

xyz_uw = read("train_UNWRAPPED.xyz",":")
mu = np.array([xyz_uw[i].info["mu"] / len(xyz_uw[i]) for i in range(len(xyz_uw))])
tmatr = np.array([[0.0,1.0,0.0],[0.0,0.0,1.0],[1.0,0.0,0.0]])
for i in range(len(mu)):
  mu[i] = np.dot(tmatr,mu[i])

mu_train = mu[:600].reshape(-1)
mu_test = mu[600:]

k_mn_data = np.dot(k_nm_train.T,mu_train)

weights = scipy.linalg.solve(kmnnm_reg,k_mn_data)
prediction = np.dot(k_nm_test,weights).reshape((-1,3))
for i in range(len(prediction)):
  prediction[i] = np.dot(tmatr.T,prediction[i]) * len(xyz_uw[i+nt])
  mu_test[i] = np.dot(tmatr.T,mu_test[i]) * len(xyz_uw[i+nt])

print(100 * np.sqrt(np.mean((prediction - mu_test)**2)) / np.sqrt(np.mean(mu_test**2)),"%")
