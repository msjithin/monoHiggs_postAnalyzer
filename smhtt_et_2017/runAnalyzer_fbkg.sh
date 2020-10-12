#!/bin/bash
set -e  # exit when any command fails
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT
unset DISPLAY

f_exe="executable_etau_fbkg"
if [ -f "$f_exe" ]; then
    echo "$f_exe exists, removing file"
    rm $f_exe
fi
./rootcom etau_analyzer_fbkg $f_exe


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
./$f_exe /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_ext1_v1_00.root DYJetsToLL_00_ext1_fbkg_test.root $nEvents 1000 2017 MC DY1JetsToLL_00
./$f_exe /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraF_05.root SingleElectron_EraF_05_fbkg_test.root $nEvents 1000 2017 DATA SingleElectron_EraF_05
