#!/bin/bash
set -e  # exit when any command fails

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


f_exe="analyze_etau_test"
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

#./$f_exe /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraE_01.root SingleElectron_EraE_01_test.root $nEvents 1000 2017 DATA SingleMuon_EraE_01 
./$f_exe  /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY1JetsToLL_04.root DY1JetsToLL_04.root -1 1000 2017 MC DY1JetsToLL_04
end=`date +%s`
runtime=$((end-start))
echo "Runtime = $runtime"
echo "Elapsed: $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
