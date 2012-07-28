#!/bin/bash
        
MPIEXEC=mpirun
NP=4

#-compile/build local executable
make

#-run the case       


mkdir -p /temp/meyerson/MHD/data
ln -s /temp/meyerson/MHD/data data
cp BOUT.inp data/BOUT.inp
cp mhd.cxx data/mhd.cxx.ref 
cp otc.grd.nc data/otc.grd.nc

echo Running with NP = $NP       
#$MPIEXEC -np $NP ./mhd
$MPIEXEC $NP ./mhd -d /temp/meyerson/MHD/data

#-check the result
#idl runidl.pro
