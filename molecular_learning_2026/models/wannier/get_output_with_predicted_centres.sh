#!/bin/bash

for el in C O;do export ee=${el};cat MU_WANNIER_${el}.out | awk '!/C|O/{if (NF==1){natom=$1}}/C|O/{if ($1==ENVIRON["ee"]){print $4*natom,$2*natom,$3*natom}}' > prediction_${el}_cartesian.txt;done
cat frames.xyz | sed "s/pos:R:3:momenta:R:3:forces:R:3/pos:R:3:n_wannier:I:1:wannier_dist:R:3/" | awk '/Lattice/{print}!/Lattice|O|H|C/{print}/H/{print $1,$2,$3,$4,0,0,0,0}/O/{print $1,$2,$3,$4,4,0,0,0}/C/{n++;if (n%2==1){nw=4}else{nw=2};print $1,$2,$3,$4,nw,0,0,0}' > orig.xyz
cat <(cat prediction_C_cartesian.txt | awk '{print "CW",$0}') <(cat prediction_O_cartesian.txt | awk '{print "OW",$0}') orig.xyz | awk '/CW/{nc++;c[nc]=$2" "$3" "$4}/OW/{no++;o[no]=$2" "$3" "$4}/Lattice/{print}!/Lattice|OW|CW/{if (NF==1){print}else if ($1=="C"){mc++;print $1,$2,$3,$4,$5,c[mc]}else if ($1=="O"){mo++;print $1,$2,$3,$4,$5,o[mo]}else{print $1,$2,$3,$4,$5,$6,$7,$8}}' > output_with_predicted_centres.xyz
