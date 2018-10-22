#!/bin/sh

read MaxIter
for i in `seq $MaxIter`
do

./bin/Release/parallel_genetic_algorithm
for num in `seq 100`
do
    if [ -e ./results/n$num ]; then
    	echo n$num exists.
    else
	mkdir ./results/n$num
	mv ./results/*.csv ./results/n$num
        echo n$num creates.
        break
    fi
done

done
