#!/bin/bash
set -e  # exit when any command fails
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

if [ -f "analyze_mutau" ]; 
then 
    echo "The file analyze_mutau exists" 
    rm analyze_mutau
fi
./rootcom mutau_analyzer analyze_mutau

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


./analyze_mutau /hdfs/store/user/jmadhusu/with_boostedtau/2017_skimmed/with_boostedtaus/mutau/DYJetsToLL_M-50_TuneCP5_ext1_v1_00/DYJetsToLL_M-50_TuneCP5_ext1_v1_00_5.root DYJetsToLL_00_test.root $nEvents 1000 2017 MC DY1JetsToLL_00
./analyze_mutau /hdfs/store/user/jmadhusu/with_boosted_taus/2017/zprimeBaryonic/Signal_Zpbaryonic2017_01_7.root Zpbaryonic2017_01_7.root $nEvents 1000 2017 MC Zpbaryonic2017_7



end=`date +%s`
runtime=$((end-start))
echo "Runtime = $runtime"
echo "Elapsed: $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
