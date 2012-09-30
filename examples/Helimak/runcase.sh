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

NOUTS=(400 200 200 400 200 300 200 100)
#tstep=(2e1 2e1 1e1 1e1 1e1 1e1 1e1) #helimak cases

# tstep=(1e6 1e5 1e5 1e7 1e5 1e7 1e5 1e5) #not drift case
# llist=(1e-5 1e-2 1e-3 1e-4 5e-3 5e-4 5e-5 1e-6) #not drift case

# tstep=(1e6 1e5 1e5 1e7 1e5 1e7 1e5 1e5) #not drift case
# llist=(1e-5 1e-2 1e-3 1e-4 5e-3 5e-4 5e-5 1e-6) #not drift case
#tstep=(1e2)
#llist=(1e-3)
#tstep=(1.8e6) #drift case
#llist=(1e-4) #drift case
#tstep=(1e6 1e6 1e6 4e5 4e5 1e4) #add rho convection
#llist=(1e-5 1e-6 1e-7 1e-4 1e-3 1e-2) #add rho convection
#NOUTS=(300 300 300 300 300 300 300 300)
NOUTS=(300 300 100 100 100 100 100 100)
#tstep=(1e7 5e7 5e7 5e7 5e7 5e7 5e7 5e7)
#tstep=(1e4 1e4 1e4 5e4 1e4 1e4 1e4 1e4) #bz only case 
#tstep=(1e0 1e0 1e0 1e0 1e0 1e0 1e0 1e0) #only ve dot grad ne term
#tstep=(1e1 1e1 1e1 1e1 1e1 1e1 1e1 1e1)
#tstep=(5e0 5e0 5e0 5e0 5e0 5e0 5e0 5e0)
tstep=(5e-1 1e0 1e0 1e-1 1e-1 1e-1 1e-1 1e-1)

llist=(1e-4 5e-3 1e-3 5e-2 1e-2 1e-3)
llist=(1e0 5e-2 5e-5)
peaklist=(1e0)
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


#rm status.log
i=0
key='haswak+lownu'
for lval in ${peaklist[@]}
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

  cp hlmk.cxx   $PWD/data_${lval}/hlmk.cxx.ref

  sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT_$key.inp > temp.inp
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
  #echo $current_dir >> status.log
  #ibrun -n $NP -o 0  ./2fluid 
  #wait
  rm -f data
done

#-check the result
#idl runidl.pro
