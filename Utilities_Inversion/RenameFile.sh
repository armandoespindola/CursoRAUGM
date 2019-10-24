#!/bin/sh
# Program Rename Files -> OBS

#Find Event Folders
folder=($(ls -l | grep shot | awk '{print $9}'));
main_fold=$(pwd)

for ifold in ${folder[@]}
do
    cd ${ifold}/DATA
    #rm *-OBS.bin
    rm *-ADJ*.bin
    
    stat=($(ls -l X*.bin | awk '{print $9}' | awk -F. '{print $1}' | awk -F- '{print $1}'))
    for istat in ${stat[@]}
    do
	echo ${ifold} ${istat}:

	cp ${istat}-VX.bin ${istat}-VX-M0.bin
	cp ${istat}-VY.bin ${istat}-VY-M0.bin
	cp ${istat}-VZ.bin ${istat}-VZ-M0.bin
    done
    cd $main_fold
    
done
