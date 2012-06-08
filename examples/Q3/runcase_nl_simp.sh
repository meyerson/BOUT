#!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    NP=256
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


llist=( 0.10 )

for lval in ${llist[@]}
do
 
  rm data_nl_s_${lval}

  mkdir data_nl_s_${lval}

  cp hlmk_nl_s+.cxx  data_nl_s_${lval}/hlmk.cxx.ref
  
  cp BOUT_nl_p.inp data_nl_s_${lval}/BOUT.inp

  #cp  data_${lval}/* data #copy the simulation data from the first part of the run
  
  sed "s/0.10_lam_n.nc/${lval}_lam_n.nc/g" BOUT_nl_p.inp > data_nl_s_${lval}/BOUT.inp
      
  ibrun -n $NP -o 0  ./hlmk_nl_s+ -d data_nl_s_${lval}
done

#-check the result
#idl runidl.pro
