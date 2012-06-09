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
data_dir='/tmp/q3'


llist=(.2 .05 .01)
#llist=('short')


rm status.log

for lval in ${llist[@]}
do
    current_dir=$data_dir/data_${lval}
    echo $current_dir
    
    rm -r $current_dir
    mkdir -p $current_dir
    
    rm -r $PWD/data_${lval}
    #mkdir $PWD/data_${lval}

    
    #rm -r data_${lval}
    
    #mkdir data_${lval}

    ltime=.001
  
    cp q3_simp.cxx   $current_dir/q3_simp.cxx.ref
    cp q3_simp.cxx   $PWD/data_${lval}/q3_simp.cxx.ref
    
    sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT.inp > temp.inp
    sed "s/TIMESTEP = .02/TIMESTEP =  ${ltime}/g" temp.inp > $current_dir/BOUT.inp

    
  #cp  data_${lval}/* data #copy the simulation data from the first part of the run
  
  #sed "s/0.10_lam_n.nc/${lval}_lam_n.nc/g" BOUT_nl_p.inp > data_nl_s_${lval}/BOUT.inp
    
    
    
    echo $MPIEXEC
    
    ln -s $current_dir $PWD/data_${lval}
    echo $current_dir >> status.log
    if [ "${lval}" = "${llist[${#llist[@]}-1]}" ]; then
	echo 'done' >> status.log
    fi

    $MPIEXEC $NP ./q3_simp -d $current_dir
    #mv post_bout.pkl $current_dir

done


#-check the result
#idl runidl.pro
