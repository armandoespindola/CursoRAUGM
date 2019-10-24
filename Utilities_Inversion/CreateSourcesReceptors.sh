#!/bin/sh
# Program creates sources and receptors files for every event and folders
# Tue Apr 16 16:18:58 +03 2019


#Find Event Folders
folder=$(echo shot{1..12})
main_fold=$(pwd)

srcline1=$(sed '1q;d' sources.par)
srcline2=$(echo "$(sed '2q;d' sources.par)" | awk -F, '{$1=$2=$3="";print $0}')

line2=$(echo $srcline2 | sed 's/ /,/g')
#echo $line2
nline=1
for ifold in ${folder}
do
    mkdir -p ${ifold}
    mkdir -p ${ifold}/{DATA,temp}

    printf "%s\n" "${srcline1}" > ${ifold}/sources.par

    line1=$(echo $(sed "${nline}q;d" src.txt) | sed 's/ /,/g')

    printf "%s,%s" $line1 $line2 >> ${ifold}/sources.par

    cp receptors.par ${ifold}/
    nline=$((nline+1))
  

done
