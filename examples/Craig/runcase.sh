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



make

rm status.log
#rm run.log

#bckey='calibrate'
bckey='dirichlet'
codekey='craig'

current_dir=$PWD
data_dir=/tmp/${codekey}

current_dir=$data_dir/data_${bckey}
echo $current_dir

#rm -r $current_dir
mkdir -p $current_dir

rm -r $PWD/data
  
  
cp ${codekey}.cxx   $current_dir/physics_code.cxx.ref
cp BOUT.inp $current_dir
 
echo $MPIEXEC
    
ln -s $current_dir $PWD/data
echo $current_dir >> status.log
    #echo $current_dir


$MPIEXEC $NP ./craig -d $current_dir
  