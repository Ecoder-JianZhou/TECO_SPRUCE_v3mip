export LD_LIBRARY_PATH=/home/jz964/miniconda3/lib/:$LD_LIBRARY_PATH 
gfortran -g -fbacktrace -Wall -fcheck=all src/dataType.f90 src/updateAndSummary.f90 src/writeOutputs2nc.f90 src/soil.f90 src/vegetation.f90 src/transfer.f90 src/driver.f90 src/mcmc.f90 src/spinup.f90  src/main.f90 -o run_teco -I/home/jz964/miniconda3/include -L/home/jz964/miniconda3/lib -lnetcdff -lnetcdf
# gfortran -g dataType.f90 driver.f90  mcmc.f90 spinup.f90 vegetation.f90 main.f90 -o run_teco
rm *.mod
./run_teco
rm run_teco
