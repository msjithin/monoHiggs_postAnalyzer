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


./rootcom smhet_2018 executable_smhtt_et
./executable_smhtt_et /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/2018/DY.root theirs_rootfile/DY.root $nEvents 1000 2018 MC DY_v1
./executable_smhtt_et /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/2018/DY1.root theirs_rootfile/DY1.root $nEvents 1000 2018 MC DY1_v1
./executable_smhtt_et /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/2018/DY2.root theirs_rootfile/DY2.root $nEvents 1000 2018 MC DY2_v1
./executable_smhtt_et /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/2018/DY3.root theirs_rootfile/DY3.root $nEvents 1000 2018 MC DY3_v1
./executable_smhtt_et /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/2018/DY4.root theirs_rootfile/DY4.root $nEvents 1000 2018 MC DY4_v1
./executable_smhtt_et /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/2018/TTTo2L2Nu.root theirs_rootfile/TTTo2L2Nu.root $nEvents 1000 2018 MC TTTo2L2Nu
./executable_smhtt_et /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/2018/TTToHadronic.root theirs_rootfile/TTToHadronic.root $nEvents 1000 2018 MC TTToHadronic
./executable_smhtt_et /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/2018/TTToSemiLeptonic.root theirs_rootfile/TTToSemiLeptonic.root $nEvents 1000 2018 MC TTToSemiLeptonic



./executable_smhtt_et /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/2018/Data.root theirs_rootfile/Data.root $nEvents 1000 2018 DATA Data
