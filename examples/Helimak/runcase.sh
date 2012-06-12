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
#cp BOUT_relax.inp BOUT.inp

#llist=(0.10 -0.10 )

current_dir=$PWD
data_dir='/tmp/hlmk'


#llist=(0.10)

NOUTS=(200 200 20 200 200 200 200 200)
tstep=(1e3 1e3 1e3 1e2 1e2 1e2 1e2 1e1)

llist=(1 1e-1 1e-2 1e-3 1e-5 1e-4)

tstep=(1e2 1e2)
llist=(1e-3 1e-4)

#rm status.log
i=0
for lval in ${llist[@]}
do
  mkdir data_${lval}
  ln -s data_${lval} data
  
  current_dir=$data_dir/data_simp_${lval}
  echo $current_dir
    
  rm -r $current_dir
  mkdir -p $current_dir
  
  rm -r $PWD/data_${lval}

  cp hlmk.cxx   $current_dir/2fluid.cxx.ref
  cp hlmk.cxx   $current_dir/hlmk.cxx.ref

  cp hlmk.cxx   $PWD/data_${lval}/hlmk.cxx.ref

  sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT.inp > temp.inp
  sed "s/NOUT = 100/NOUT = ${NOUTS[$i]}/g" temp.inp > temp2.inp
  sed "s/TIMESTEP = 5e2/TIMESTEP =  ${tstep[$i]}/g" temp2.inp > $current_dir/BOUT.inp

#  cp hlmk.cxx data/2fluid.cxx.ref
 # cp hlmk.cxx data/hlmk.cxx.ref #copy the source code for easy reference
  #cp  data_${lval}/* data #copy the simulation data from the first part of the run
  
  #sed "s/..\/Helimak\/Ln_lam\/Helimak_32x32_0.10_lam_n.nc/..\/Helimak\/Ln_lam\/Helimak_32x32_${lval}_lam_n.nc/g" BOUT.inp > data/BOUT.inp
  #sed "s/Helimak_4x32_0.10_lam_n.nc/Helimak_4x32_${lval}_lam_n.nc/g" BOUT.inp >  $current_dir/BOUT.inp
  #if [ $zval -lt 128 ]
      #then
      # reduce time-step. At large times these cases produce noise
      #sed "s/TIMESTEP = 5e3/TIMESTEP = 1e3/g" data/BOUT.inp > data/tmp
      #mv -f data/tmp data/BOUT.inp
  #fi
  #cp BOUT.inp $current_dir/BOUT.inp
  echo "$((i++))"  
  
  $MPIEXEC $NP ./hlmk -d $current_dir

  ln -s $current_dir $PWD/data_${lval}
  echo $current_dir >> status.log
  #ibrun -n $NP -o 0  ./2fluid 
  #wait
  rm -f data
done

#-check the result
#idl runidl.pro
