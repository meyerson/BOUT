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


llist=('hd')

for lval in ${llist[@]}
do
 
  #rm -r data_${lval}

  #mkdir data_re_${lval}
  
  rm data_re_${lval}/*

  cp data_${lval}/BOUT.restart_0100* data_re_${lval}/

  cd data_re_${lval}

  rename s/t_0100/t/g *.nc

  cd ..

  cp q3_simp_re.cxx  data_re_${lval}/q3_simp.cxx.ref
  
  cp BOUT_re_hd.inp data_re_${lval}/BOUT.inp

  #cp  data_${lval}/* data #copy the simulation data from the first part of the run
  
  #sed "s/0.10_lam_n.nc/${lval}_lam_n.nc/g" BOUT_nl_p.inp > data_nl_s_${lval}/BOUT.inp
      
  $MPIEXEC $NP ./q3_simp_re re -d data_re_${lval}
done

#-check the result
#idl runidl.pro
