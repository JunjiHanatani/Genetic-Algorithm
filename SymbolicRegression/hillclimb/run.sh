#!/bin/sh

echo "Input the number of the runs."
read MaxIter
for i in `seq $MaxIter`
do

./bin/Release/hillclimb
for num in `seq 100`
do
    folder=`expr $num - 1`
    if [ -e ./results/n$folder ]; then
    	echo n$folder exists.
    else
	mkdir ./results/n$folder
	mv ./results/*.csv ./results/n$folder
        echo n$folder creates.
        break
    fi
done

done
