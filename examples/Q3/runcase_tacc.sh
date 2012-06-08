
#!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    NP=16
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
echo $current_dir

#data_dir=${current_dir/"/share/home/01523"/"/scratch/01523"}
data_dir=`echo $PWD|sed "s/share\/home/scratch/g"`


llist=(1 .5 .1 .05)

for lval in ${llist[@]}
do
    current_dir=$data_dir/data_${lval}
    echo $current_dir

    rm -r $current_dir
    mkdir $current_dir

  
  #ltime=$(echo "scale=5; ${lval}*.02" | bc -l)
    ltime=.01
    
    echo $ltime
    
    cp q3_simp.cxx $current_dir/q3_simp.cxx.ref
    
    sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT.inp > temp.inp
    sed "s/TIMESTEP = 1/TIMESTEP =  ${ltime}/g" temp.inp > $current_dir/BOUT.inp

    ibrun -n $NP -o 0  ./q3_simp -d $current_dir

done

#-check the result
#idl runidl.pro
