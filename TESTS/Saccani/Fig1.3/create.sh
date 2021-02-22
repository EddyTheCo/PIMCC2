#!/bin/bash
for landa in   0.5 ;
do

for beta in   1.0; 
do
for bound in  periodic;
do
                for tao in 0.0078125 0.00390625;
                do
                        for pot in free;
                        do
                                for opt in false ;
                                do
                                      for Npar in  1 ;
                                      do
                        for ensem  in  GrandCanonical ;
do
for mu in  37 38 39 40 41 42;
do
for Bpot in 5.0;
do
for eta in 1;
do
                                TimeSlices=$(echo $beta  $tao= | awk '{ print $1/$2 }')
#                               MBar=$(echo  $TimeSlices= | awk '{ print $1/2 }')
MBar=32;
                                 mkdir Beta.$beta.N.$Npar.NT.$TimeSlices.Sacc.mu.$mu.CB.$Bpot.eta.$eta
                        cd Beta.$beta.N.$Npar.NT.$TimeSlices.Sacc.mu.$mu.CB.$Bpot.eta.$eta

cat > input  <<EOF
$Npar              #OfParticles, doesnt do nothing if restart
$tao         #Tao
$beta              #beta
$landa              #landa
100000            #OfBlocks
2               #dimension
11.855 10.267  #L d numbers
100           #NumberOfsweepsInABlock
$mu            #mu
start           #start or restart
$ensem       # or Canonical
0 0       #K d numbers from (k1*x1^2+..)/2
$eta               #eta factor
periodic        #free or periodic
$MBar           #MBar
$pot           #free or harmonic or infinteWell
LiebLiniger           #LiebLiniger or free or Dipolar or Softcore or Bonin
$Bpot              #CLieb or CDip or V_0Softcore or SigmaBonin
false               #FourOrderPotential true or false
4                   #of Sampling
0               #Warmup
EOF

cd ../
done
done
done
done
done
done
done
done
done
done
done

