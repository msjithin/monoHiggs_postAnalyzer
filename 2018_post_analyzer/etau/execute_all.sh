
set -e

if [ "$(ls -A "rootFiles")" ]; then
    rm rootFiles/*
fi
cp /nfs_scratch/jmadhusu/CMSSW_10_2_10/src/2018_skim/etau/analyzer/output/*.root rootFiles/

sh hadd_files.sh

cd plotting_script

sh _postAnalyzer_mutau.sh
