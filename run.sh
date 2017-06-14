#!/bin/bash

name=$(sed -n "6 p" poisson.in)
#echo -n "Write data to folder > "
#read name
#./poisson -ksp_type gmres -log_view #-ksp_monitor_max #-mat_view ::ascii_matlab
#./poisson -pc_type jacobi -ksp_type richardson -ksp_richardson_scale 0.66666666 -log_view -ksp_monitor_max #-mat_view ::ascii_matlab
./poisson -pc_type jacobi -ksp_type gmres -log_view -ksp_monitor_max #-mat_view ::ascii_matlab
mkdir -p data/$name
#mv -v *.out data/$name/
cp -r *.in data/$name/
mv -v *.err data/$name/
mv -v *.dat data/$name/


