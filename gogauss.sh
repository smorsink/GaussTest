#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times

#base="/home/kitung"
base="/Users/sharon/code"
exe_dir="$base/Marginalize-master"

#make gauss1d

data="$exe_dir/dataGausstest"

# ferret version 1e2 points
#gauss1d -i "$data/data1d1e2.txt" -o "out_mean1de2.txt" -O "out_sig1de2.txt" -n 100 -b 400 -p 100 -c 100000000 -f 0.1 -a "$exe_dir/1d1e2/data"

# ferret version 1e4 points -- FINAL VERSION
#gauss1d -i "$data/data1d1e4.txt" -o "out_mean1de4.txt" -O "out_sig1de4.txt" -n 10000 -b 100000 -p 8000 -c 100000000 -f 0.01 -a "$exe_dir/1d1e4/data"

# pure metropolis-hastings
#gauss1d -i "$data/data1d1e2.txt" -o "out1d1e2/out_mean1de2.txt" -O "out1d1e2/out_sig1de2.txt" -n 100 -b 1600 -c 100000 -f 0.1 -p 1600

#Good values for 10^4 data set
#gauss1d -i "$data/data1d1e4.dat" -o "out_mean1de4.txt" -O "out_sig1de4.txt" -n 10000 -b 100000 -c 1000000 -f 0.01 -p 100000

#make gauss2d -i "$data/data2d1e4.dat" -o "out_mean2de4.txt" -O "out_sig2de4.txt" -n 10000 -b 2000 -c 10000000 -f 0.2

#gauss2d 

make gauss2d


#1e2 Version
#gauss2d -i "$data/data2d1e2.txt" -n 100 -b 1000 -p 1000 -f 1.0 -c 10000

#1e4 Version
gauss2d -i "$data/data2d1e4.dat" -n 10000 -b 40 -p 100 -f 1.0 -c 100




times

