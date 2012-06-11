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
make

#-run the case       
echo Running with NP = $NP       

rm -rf data*
cp BOUT_relax.inp BOUT.inp

#llist=(0.10 -0.10 )
llist=(0.10)

for lval in ${llist[@]}
do
  mkdir data_${lval}
  ln -s data_${lval} data
  
  cp hlmk.cxx data/2fluid.cxx.ref
  cp hlmk.cxx data/hlmk.cxx.ref #copy the source code for easy reference
  #cp  data_${lval}/* data #copy the simulation data from the first part of the run
  
  #sed "s/..\/Helimak\/Ln_lam\/Helimak_32x32_0.10_lam_n.nc/..\/Helimak\/Ln_lam\/Helimak_32x32_${lval}_lam_n.nc/g" BOUT.inp > data/BOUT.inp
  sed "s/Helimak_4x32_0.10_lam_n.nc/Helimak_4x32_${lval}_lam_n.nc/g" BOUT.inp > data/BOUT.inp
  #if [ $zval -lt 128 ]
      #then
      # reduce time-step. At large times these cases produce noise
      #sed "s/TIMESTEP = 5e3/TIMESTEP = 1e3/g" data/BOUT.inp > data/tmp
      #mv -f data/tmp data/BOUT.inp
  #fi
      
  $MPIEXEC $NP ./hlmk
  #ibrun -n $NP -o 0  ./2fluid 
  #wait
  rm -f data
done

#-check the result
#idl runidl.pro
