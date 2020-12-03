#!/bin/bash
set -e  # exit when any command fails

nEvents=-1
sample='dy'
plottingOn=0
while getopts n:s:p option
do
    case "${option}"
	in
	n) nEvents=${OPTARG};;
	p) plottingOn=1 ;;
	s) sample=${OPTARG};;
esac
done


./rootcom smhet_2017_fbkg executable_smhtt_et_fbkg
./executable_smhtt_et_fbkg ../theirs_rootfile/fromCecile_DY.root theirs_DY_fbkg.root $nEvents 1000 2017 MC DY_v1
#./executable_smhtt_et_fbkg theirs_rootfile/fromCecile_Data.root theirs_rootfile/Data_fbkg.root $nEvents 1000 2017 DATA Data
