#!/bin/bash
set -e  # exit when any command fails

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

./rootcom skimm_tt_2018 analyze_tautau_updated

outFile="outSkimmed_mt.root"
nEvents=100000
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

outFile="outSkimmed_tt_dy.root"
echo "dy sample analysis....."
./analyze_tautau_updated /hdfs/store/user/jmadhusu/MC2018_Autumn18_monoHiggs_09Jun2020/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_DY1JetsToLL/200610_213927/0000/ $outFile $nEvents 1000 2018_test MC

echo "out put written to $outFile"
outFile="outSkimmed_tt_data.root"
echo "data sample analysis....."
./analyze_tautau_updated /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauA/200609_122655/0000/ $outFile $nEvents 1000 2018_test DATA
echo "out put written to $outFile"


