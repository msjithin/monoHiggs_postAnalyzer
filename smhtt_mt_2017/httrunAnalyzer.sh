#!/bin/bash
set -e  # exit when any command fails

nEvents=-1
sample='dy'
plottingOn=0
condorSubmit=0
while getopts n:s:pc option
do
    case "${option}"
	in
	n) nEvents=${OPTARG};;
	p) plottingOn=1 ;;
	s) sample=${OPTARG};;
	c) condorSubmit=1 ;;
esac
done

if (( condorSubmit == 1 ))
then 
    echo "condorSubmit turned ON"
    outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
    mkdir $outDir 
    
    ###########################   MC  #########################
    
    ./rootcom smhet_2017 executable_smhtt_mt
    
    ./MakeCondorFiles.csh executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY.root DY.root -1 1000 2017 MC DY_v1 $outDir
    ./MakeCondorFiles.csh executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY1.root DY1.root -1 1000 2017 MC DY1_v1 $outDir
    ./MakeCondorFiles.csh executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY2.root DY2.root -1 1000 2017 MC DY2_v1 $outDir
    ./MakeCondorFiles.csh executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY3.root DY3.root -1 1000 2017 MC DY3_v1 $outDir
    ./MakeCondorFiles.csh executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY4.root DY4.root -1 1000 2017 MC DY4_v1 $outDir
    ./MakeCondorFiles.csh executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/TTTo2L2Nu.root TTTo2L2Nu.root -1 1000 2017 MC TTTo2L2Nu $outDir
    ./MakeCondorFiles.csh executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/TTToHadronic.root TTToHadronic.root -1 1000 2017 MC TTToHadronic $outDir
    ./MakeCondorFiles.csh executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/TTToSemiLeptonic.root TTToSemiLeptonic.root -1 1000 2017 MC TTToSemiLeptonic $outDir
    ./MakeCondorFiles.csh executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/Data.root Data.root -1 1000 2017 DATA Data $outDir
    
else
    echo "condorSubmit turned OFF"
    ./rootcom smhet_2017 executable_smhtt_mt
    ./executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY.root theirs_rootfile/DY.root $nEvents 1000 2017 MC DY_v1
    ./executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY1.root theirs_rootfile/DY1.root $nEvents 1000 2017 MC DY1_v1
    ./executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY2.root theirs_rootfile/DY2.root $nEvents 1000 2017 MC DY2_v1
    ./executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY3.root theirs_rootfile/DY3.root $nEvents 1000 2017 MC DY3_v1
    ./executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/DY4.root theirs_rootfile/DY4.root $nEvents 1000 2017 MC DY4_v1
    ./executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/TTTo2L2Nu.root theirs_rootfile/TTTo2L2Nu.root $nEvents 1000 2017 MC TTTo2L2Nu
    ./executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/TTToHadronic.root theirs_rootfile/TTToHadronic.root $nEvents 1000 2017 MC TTToHadronic
    ./executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/TTToSemiLeptonic.root theirs_rootfile/TTToSemiLeptonic.root $nEvents 1000 2017 MC TTToSemiLeptonic
    ./executable_smhtt_mt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/mutau2017/Data.root theirs_rootfile/Data.root $nEvents 1000 2017 DATA Data
    
fi


