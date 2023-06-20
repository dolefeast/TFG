#! /bin/bash

mkdir logfiles_phase2
for dir in * ;do 
	echo $dir
	cp $dir/log*2nd* logfiles_phase2
done
