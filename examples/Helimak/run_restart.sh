#!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    NP=4
    MPIEXEC="mpirun -np"

fi  
# Usage: scriptname -options
# Note: dash (-) necessary

while getopts ":n:np" Option
do
  case $Option in
    n ) MPIEXEC="mpirun -np";NP=$OPTARG;;
    * ) ;;   # DEFAULT
  esac
done

#-compile/build local executable
#make

#-run the case       
echo Running with NP = $NP       

#rm -rf data

#ilist=(0 1 2 3 4 5 6 7 8 9)

#mkdir data

#for ival in ${ilist[@]}
#do
#  mkdir ./H_Drift/data_${ival}
#  ln -s ./H_Drift/data_${ival} data
  #sed "s/Zeff = 128.0/Zeff = ${zval}.0/g" BOUT.inp > data/BOUT.inp

cp BOUT.inp data/BOUT.inp

#  sed "s/10x64_grid.nc/10x64_grid${ival}.nc/g" BOUT_drift.inp > data/BOUT.inp
  #if [ $ival -lt 128 ]
  #    then
      # reduce time-step. At large times these cases produce noise
      #sed "s/TIMESTEP = 5e3/TIMESTEP = 1e3/g" data/BOUT.inp > data/tmp
      
   #   mv -f data/tmp data/BOUT.inp
  #fi
      
  #$MPIEXEC $NP ./2fluid
  ibrun -n $NP -o 0  ./2fluid re
  wait
  #rm -rf data
done

#-check the result
#idl runidl.pro
