#!/bin/bash

name=$(sed -n "6 p" poisson.in)
#echo -n "Write data to folder > "
#read name
#./poisson -pc_type ilu -ksp_type gmres -ksp_monitor_max -log_view #-mat_view ::ascii_matlab
#./poisson -pc_type ilu -ksp_type gmres -info #-mat_view ::ascii_matlab
#./poisson -pc_type ilu -ksp_type gmres -ksp_monitor_max #-mat_view ::ascii_matlab
#./poisson -pc_type none -ksp_type gmres -ksp_monitor_max #-mat_view ::ascii_matlab
#./poisson -pc_type jacobi -ksp_type richardson -ksp_richardson_scale 0.66666666 -ksp_monitor_max #-mat_view ::ascii_matlab
#./poisson -pc_type jacobi -pc_jacobi_type diagonal -ksp_type richardson -ksp_richardson_scale 1.0 -ksp_monitor_max #-mat_view ::ascii_matlab
#./poisson -pc_type jacobi -ksp_type gmres -ksp_monitor_max #-mat_view ::ascii_matlab
#./poisson -ksp_type bicg -ksp_monitor_max #-mat_view ::ascii_matlab
#./poisson -ksp_type bcgs -ksp_monitor_max #-mat_view ::ascii_matlab
#./poisson -pc_type ilu -ksp_type gmres -ksp_monitor_draw #-mat_view ::ascii_matlab


#mpiexec -n 2 ./poisson -pc_type bjacobi -ksp_type gmres -ksp_view -ksp_monitor -log_view #-mat_view ::ascii_matlab
#mpiexec -n 4 ./poisson -ksp_type gmres -ksp_view -ksp_monitor -log_view #-mat_view ::ascii_matlab
mpiexec -n 1 ./poisson -ksp_type gmres -pc_asm_blocks 32 -ksp_view #-mat_view ::ascii_matlab
#mpiexec -n 8 ./poisson -pc_type bjacobi -ksp_type gmres -ksp_view -log_view #-mat_view ::ascii_matlab
#mpiexec -n 8 ./poisson -pc_type jacobi -ksp_type richardson -ksp_richardson_scale 0.66666666 -ksp_view -ksp_monitor -log_view #-mat_view ::ascii_matlab

#mpiexec -n 2 ./poisson -log_view #-mat_view ::ascii_matlab
mkdir -p data/$name
mv -v *.out data/$name/
cp -r *.in data/$name/
cp -r *.sh data/$name/
mv -v *.err data/$name/
mv -v *.dat data/$name/


