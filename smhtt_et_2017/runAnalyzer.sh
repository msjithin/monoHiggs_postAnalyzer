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
./executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_ext1_v1_00.root DYJetsToLL_00_ext1_test.root $nEvents 1000 2017 MC DY1JetsToLL_00
#./executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraF_05.root SingleElectron_EraF_05.root $nEvents 1000 2017 DATA SingleElectron_EraF_05
