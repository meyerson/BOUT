#!/bin/bash

    NP=4
    MPIEXEC="mpirun -np"

echo data$1
if [ $# = "0" ]; then
working_dir = data
fi

mkdir data$1

cp BOUT.inp data2/BOUT.inp


#ibrun -n $NP -o 0  ./2fluid
wait
  #rm -rf data
#done

#-check the result
#idl runidl.pro
