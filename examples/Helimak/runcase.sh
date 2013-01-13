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


#NOUTS=(300 300 300 300 300 100 100 100)
NOUTS=(100 200 100 200 200 200)

#tstep=(1e5 1e-1 1e1 1e-3 1e-2 1e-3) #bz only steps

tstep=(1e2 1e-1 1e1 1e-1 1e-1 5e-1) #bz_1_10 steps

#llist=(1e-4 5e-3 1e-3 5e-2 1e-2 1e-1) #bz and bz11


llist=(1e-3 5e-2 1e-2 5e-1 1e-1 1e0) #bz1_10

NOUTS=(200 200)
tstep=(1e0 1e0)
llist=(1e-1 1e0)

# NOUTS=(100)
# tstep=(1e0)
# llist=(1e0)

#works with TOL 1e-8
#tstep=(1e-3)
#llist=(1e-1)

#works with TOL 1e-8
#tstep=(1e1 1e-3)
#llist=(1e-3 1e-1)

#tstep=(3e2)
#llist=(1e-4)


#peaklist=(1e0)
#tstep=(3e1 3e1 3e1 3e1 3e1 3e1 3e1 3e1)
#tstep=(1e1 1e1 1e1 1e1 1e1 1e1 1e1 1e1)
#tstep=(1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1)
#tstep=(1e2 1e2 1e2 1e2 1e2 1e2 1e2 1e2)
#tstep=(1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1)
#tstep=(3e-2 3e-2 3e-2 3e-2 3e-2 3e-2 3e-2 3e-2)
restep=(2e-4 1e-4 5e-4 4e-4 3e-4 4e-4)
detlist=(1e-3 5e-4 1e-4 5e-3 4e-4 3e-4 2e-4)

#tstep=(1e0 1e0)
#detlist=(1e-3 5e-4)

tstep2=(1e-2 1e-2 1e-2 1e-2)
shiftlist=(1.7e-3 5.7e-4 7e-4)

#tstep2=(1e-3 1e-3 1e-3 1e-3 5e-4 5e-4)
#shiftlist=(1e-2 5e-3 1e-3 3.7e-4 5.7e-5 1e-4)



i=0
#key='haswak+lownu+boost+bz1_10'
#key='haswak+lownu+4fld+bz1_10_euler'
key='haswak'

rm status_${key}.log

for lval in ${llist[@]}
do
  mkdir data_${lval}
  ln -s data_${lval} data
  
  current_dir=$data_dir/data_bz_${key}_${lval}
  echo $current_dir
    
  #rm -r $current_dir
  mkdir -p $current_dir
  
  rm -r $PWD/data_${lval}

  cp hlmk.cxx   $current_dir/2fluid.cxx.ref
  cp hlmk.cxx   $current_dir/hlmk.cxx.ref
  cp hlmk.cxx   $current_dir/physics_code.cxx.ref

  cp hlmk.cxx   $PWD/data_${lval}/hlmk.cxx.ref
  #cp $data_dir/data_bz_${key}_${lval}/*restart* $current_dir

  sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT_${key}.inp > temp.inp
  #sed "s/ZMAX = 1/ZMAX = 1/g" BOUT.inp > temp.inp
  sed "s/NOUT = 100/NOUT = ${NOUTS[$i]}/g" temp.inp > temp2.inp
  #sed "s/NOUT = 100/NOUT = 100/g" temp.inp > temp2.inp
#sed "s/TIMESTEP = 5e2/TIMESTEP =  ${tstep[$i]}/g" temp2.inp > $current_dir/BOUT.inp
  sed "s/TIMESTEP = 5e2/TIMESTEP =  ${tstep[$i]}/g" temp2.inp > $current_dir/BOUT.inp
  #sed "s/TIMESTEP = 5e2/TIMESTEP =  5e-2/g" temp2.inp > $current_dir/BOUT.inp
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
  echo $current_dir >> status_${key}.log
  #ibrun -n $NP -o 0  ./2fluid 
  #wait
  rm -f data
done

#-check the resulthl
#idl runidl.pro
