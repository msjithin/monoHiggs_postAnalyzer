#!/bin/bash
set -e  # exit when any command fails

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

if [ -f "analyze_tautau" ]; 
then 
    echo "The file analyze_tautau exists" 
    rm analyze_tautau
fi
./rootcom tautau_analyzer analyze_tautau

outFile="study_mutau_110k.root"
start=`date +%s`
nEvents=1000
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
./analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau/DYJetsToLL_M-50_TuneCP5_v1_00.root DYJetsToLL_00_test.root $nEvents 1000 2017 MC DY1JetsToLL_00

echo "wjets sample analysis....."
./analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau/WJetsToLNu_TuneCP5_00.root WJetsToLNu_00_test.root $nEvents 1000 2017 MC WJetsToLNu_00

echo "data sample analysis....."
./analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau/Tau_EraC_01.root Tau_EraC_01_test.root $nEvents 1000 2017 DATA Tau_EraC_01

end=`date +%s`
runtime=$((end-start))
echo "Runtime = $runtime"
echo "Elapsed: $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
