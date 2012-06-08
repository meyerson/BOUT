#!/bin/bash

llist=(1 .5 .1 .05)
for lval in ${llist[@]}
do
    echo ${lval}
    echo ${llist[${#llist[@]}-1]}
    
    if [ "${lval}" = "${llist[${#llist[@]}-1]}" ]; then
	echo 'last one'
    fi
done