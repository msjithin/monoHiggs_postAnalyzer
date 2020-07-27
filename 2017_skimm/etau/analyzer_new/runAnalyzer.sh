#!/bin/bash
set -e  # exit when any command fails

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


f_exe="analyze_etau"
if [ -f "$f_exe" ]; then
    echo "$f_exe exists, removing file"
    rm $f_exe
fi

./rootcom etau_analyzer $f_exe


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
./$f_exe /hdfs/store/user/jmadhusu/2017_skimmed/etau/DYJetsToLL_M-50_TuneCP5_v1_00.root DYJetsToLL_00_test.root $nEvents 1000 2017 MC DY1JetsToLL_00
./$f_exe /hdfs/store/user/jmadhusu/2017_skimmed/etau/ee_DY1JetsToLL_M-50_TuneCP5_01.root ee_DY1JetsToLL_M-50_TuneCP5_01_test.root $nEvents 1000 2017 MC DY1JetsToLL_01
echo "wjets sample analysis....."
./$f_exe /hdfs/store/user/jmadhusu/2017_skimmed/etau/WJetsToLNu_TuneCP5_00.root WJetsToLNu_00_test.root $nEvents 1000 2017 MC WJetsToLNu_00

echo "data sample analysis....."
./$f_exe /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraB_00.root SingleElectron_EraB_00_test.root $nEvents 1000 2017 DATA SingleMuon_EraF_03

./$f_exe /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraE_01.root SingleElectron_EraE_01_test.root $nEvents 1000 2017 DATA SingleMuon_EraE_01 

end=`date +%s`
runtime=$((end-start))
echo "Runtime = $runtime"
echo "Elapsed: $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
