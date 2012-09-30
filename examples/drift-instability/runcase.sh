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

current_dir=$PWD
data_dir='/tmp/2fluid'

#full list
llist=(1e0 5e-1 5e-2 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7)

#short list
#llist=(1e-3 1e-4 1e-5 1e-6 1e-7)
#llist=(8e-4)
#1e-4  1e-2
#llist=(1.08e-5 2e-5)
#llist=(1e-2 1e-3 1e-1)
#llist=(1e-5)
#5e-6 1e-5 5e-5 6e-5 8e-5 1e-4 2e-4 5e-4 1e-3 2e-3 1e-2)
#tstep=(1e2 1e2 1e3)
#tstep=(1e3 2e3)

NOUTS=(100 200 500 500 200 200 200 100)
#tstep=(1e1 1e2 1e2 2e2 1e2 1e2 1e2 1e2) #large b_phi case

NOUTS=(100 100 100 100 100 200 200 100)
#tstep=(1e2 1e2 1e2 5e1 5e1 1e2 1e2 1e2) #small b_phi case
#NOUTS=(10 10 20 20 200 200 200 100)
#tstep=(2e1 2e1 1e1 1e1 1e2 1e2 1e2 1e2) #lar
#llist=(5e-5 4e-4 2e-4 2e-5 8e-3 1e-3 5e-3)

#tstep=(1e2 5e1 1e1 2e1 5e1)
tstep=(1e1 1e1 1e1 1e1 1e1 1e2 1e2)
llist=(1e-3 5e-4 1e-4 5e-5 1e-5)
#llist=(8e-3 8e-4 8e-5)
#tstep=(1e1)
#llist=(1e-5)

make

rm status.log
#rm run.log
i=0
key='AA36_grid11'
for lval in ${llist[@]}
do
    current_dir=$data_dir/data_${key}_${lval}
    echo $current_dir
    
    rm -r $current_dir
    mkdir -p $current_dir
    
    rm -r $PWD/data_${lval}
    #mkdir $PWD/data_${lval}

    
    #rm -r data_${lval}
    
    #mkdir data_${lval}
    #echo "scale=10; 4*$lval*100000.0" | bc -l > ltime
    
    
    #for slower growing domains we need to up the timestep
  
    cp 2fluid.cxx   $current_dir/2fluid.cxx.ref
    cp 2fluid.cxx   $PWD/data_${lval}/2fluid.cxx.ref
    

    sed "s/ZMAX = .0001/ZMAX = ${lval}/g" BOUT.inp > temp.inp
    sed "s/NOUT = 200/NOUT = ${NOUTS[$i]}/g" temp.inp > temp2.inp
    #sed "s/NOUT = 200/NOUT = 100/g" temp.inp > temp2.inp

    #sed "s/TIMESTEP = 5e3/TIMESTEP =  1e2/g" temp2.inp > $current_dir/BOUT.inp
    sed "s/TIMESTEP = 5e3/TIMESTEP =  ${tstep[$i]}/g" temp2.inp > $current_dir/BOUT.inp
    
    echo "$((i++))"
    
    
  #cp  data_${lval}/* data #copy the simulation data from the first part of the run
  
  #sed "s/0.10_lam_n.nc/${lval}_lam_n.nc/g" BOUT_nl_p.inp > data_nl_s_${lval}/BOUT.inp
    
    
    
    echo $MPIEXEC
    
    ln -s $current_dir $PWD/data_${lval}
    echo $current_dir >> status.log
    #echo $current_dir


    $MPIEXEC $NP ./2fluid -d $current_dir 
    #>> run.log
    #mv post_bout.pkl $current_dir
    

    

    
done

#-check the result
#idl runidl.pro
