#!/bin/bash
set -e  # exit when any command fails

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

./rootcom skimm_mt_2017 analyze_mutau_updated

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

outFile="outSkimmed_mt_dy.root"
echo "dy sample analysis....."
./analyze_mutau_updated /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_DY1JetsToLL_M-50_TuneCP5/200609_002516/0000/ $outFile $nEvents 1000 2018_test MC
outFile="outSkimmed_et_dy2.root"
./analyze_mutau_updated /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_DY2JetsToLL_M-50_TuneCP5/200609_002539/0000/ $outFile $nEvents 1000 2018_test MC
echo "out put written to $outFile"
outFile="outSkimmed_mt_data.root"
echo "data sample analysis....."
./analyze_mutau_updated /hdfs/store/user/jmadhusu/data2017_31Mar2018_09Jun2020/SingleMuon/crab_job_SingleMuon_EraB/200609_001433/0000/ $outFile $nEvents 1000 2018_test DATA
echo "out put written to $outFile"


