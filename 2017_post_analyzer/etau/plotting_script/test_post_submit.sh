#!/bin/bash 
set -e 

if [ -f "f_mutau_initial.root" ]; then
    echo "deleting existing f_mutau_initial.root file ....."
    rm f_mutau_initial.root
fi
if [ "$(ls -A files_initial)" ]; then
    echo "deleting existing files in directory files_initital ....."
    rm files_initial/*.root
fi

./Make.sh _postAnalyzer_mutau.C 
echo "making files_initial......."
./_postAnalyzer_mutau.exe ../files_nominal/DY1JetsToLL_final.root files_initial/DY1JetsToLL_final.root ZTT1jet ZTT1jet 0
./_postAnalyzer_mutau.exe ../files_nominal/DY2JetsToLL_final.root files_initial/DY2JetsToLL_final.root ZTT2jet ZTT2jet 0
./_postAnalyzer_mutau.exe ../files_nominal/DY3JetsToLL_final.root files_initial/DY3JetsToLL_final.root ZTT3jet ZTT3jet 0
./_postAnalyzer_mutau.exe ../files_nominal/DY4JetsToLL_final.root files_initial/DY4JetsToLL_final.root ZTT4jet ZTT4jet 0
./_postAnalyzer_mutau.exe ../files_nominal/DYJetsToLL_final.root files_initial/DYJetsToLL_final.root ZTTinc  ZTTinc 0
#./_postAnalyzer_mutau.exe ../files_nominal/DYJetsToLL_M10to50_final.root files_initial/DYJetsToLL_M10to50_final.root DYJetsToLL_M10to50_final DYJetsToLL_M10to50_final 0
./_postAnalyzer_mutau.exe ../files_nominal/SingleMuonA_final.root files_initial/SingleMuonA_final.root data_obs data_obs 0
./_postAnalyzer_mutau.exe ../files_nominal/SingleMuonB_final.root files_initial/SingleMuonB_final.root data_obs data_obs 0
./_postAnalyzer_mutau.exe ../files_nominal/SingleMuonC_final.root files_initial/SingleMuonC_final.root data_obs data_obs 0
./_postAnalyzer_mutau.exe ../files_nominal/SingleMuonD_PromptReco_1_final.root files_initial/SingleMuonD_PromptReco_1_final.root data_obs data_obs 0
./_postAnalyzer_mutau.exe ../files_nominal/TTTo2L2Nu_powheg_final.root files_initial/TTTo2L2Nu_powheg_final.root TTTo2L2Nu TTTo2L2Nu 0
./_postAnalyzer_mutau.exe ../files_nominal/TTToHadronic_powheg_final.root files_initial/TTToHadronic_powheg_final.root TTToHadronic TTToHadronic 0
./_postAnalyzer_mutau.exe ../files_nominal/TTToSemiLeptonic_powheg_final.root files_initial/TTToSemiLeptonic_powheg_final.root TTToSemiLeptonic TTToSemiLeptonic 0 
hadd -f f_mutau_initial.root files_initial/*.root 
echo "*************** root file made ***************" 
#sh do_plots_mutau.sh 
echo "*************** plots made ***************" 
