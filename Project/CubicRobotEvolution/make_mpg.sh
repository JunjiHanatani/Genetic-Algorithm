#!/bin/sh

for num in `seq 100`
do
    folder=`expr $num - 1`
    if [ -e ./png/n$folder ]; then
    	echo n$folder exists.
    else
	mkdir ./png/n$folder
	mv ./png/*.png ./png/n$folder
        echo n$folder creates.
	ffmpeg -framerate 120 -i ./png/n$folder/1%04d.png -vcodec libx264 -pix_fmt yuv420p -r 60 ./out_n$folder.mp4
        break
    fi
done


