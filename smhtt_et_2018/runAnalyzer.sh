#!/bin/bash
set -e  # exit when any command fails


./rootcom etau_analyzer executable_etau


outFile="study_mutau_110k.root"
start=`date +%s`
nEvents=10000
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

echo "dy sample analysis....."
./executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DYJetsToLL_00.root DYJetsToLL_00_test.root $nEvents 1000 2018 MC DY1JetsToLL_00
./executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_00.root EGamma2018A_00_test.root $nEvents 1000 2018 DATA SingleElectron_EraF_05

