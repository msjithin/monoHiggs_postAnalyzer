#!/bin/bash
set -e 
extn="*.err"
fDir="Out_MC_08-05-2020_17-48"
while getopts "eold:" option;
do
    case $option in
        e) extn="*.err";;
	o) extn="*.out";;
	l) extn="*.log";;
	d) fDir=${OPTARG};;
    esac
done

echo "checking $extn in $fDir  ..."
for f in `ls ${fDir}`; do
    if [[ $f == $extn ]]
    then
	#echo "Processing ${f} "
	if [[ -s $f ]]
	then
	     echo "$f has some data."
	fi
    fi
done

