#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times

base="/home/kitung"
#base="/Users/sharon/code"
exe_dir="$base/GaussTest"

make gauss1d

data="$exe_dir/dataGausstest"


./loglike1d -i "$data/data1d1e2.txt" -n 100 -m 0.0 -s 1.0 -b 400

#./gauss1d -i "$data/data1d1e2.txt" -o "out_mean1de2.txt" -O "out_sig1de2.txt" -n 100 -b 400 -c 1000000 -f 0.4

#Good values for 10^4 data set
#gauss1d -i "$data/data1d1e4.dat" -o "out_mean1de4.txt" -O "out_sig1de4.txt" -n 10000 -b 2000 -c 10000000 -f 0.2

times
