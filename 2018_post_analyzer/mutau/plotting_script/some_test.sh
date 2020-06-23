
#!/bin/bash 
set -e 
if [ -f "f_mutau_initial.root" ]; then
    echo "deleting existing f_mutau_initial.root file ....."
    rm f_mutau_initial.root
fi
if [ "$(ls -A files_nominal)" ]; then
    echo "deleting existing files in directory files_nominal ....."
    rm files_nominal/*.root
fi


./Make.sh _postAnalyzer_mutau.C 

./_postAnalyzer_mutau.exe ../files_initial/DY1JetsToLL_final.root files_nominal/DY1JetsToLL_final.root ZTT1 ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DY2JetsToLL_final.root files_nominal/DY2JetsToLL_final.root ZTT2 ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DY3JetsToLL_final.root files_nominal/DY3JetsToLL_final.root ZTT3 ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DY4JetsToLL_final.root files_nominal/DY4JetsToLL_final.root ZTT4 ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DYJetsToLL_final.root files_nominal/DYJetsToLL_final.root ZTTinc ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DYJetsToLL_M10to50_final.root files_nominal/DYJetsToLL_M10to50_final.root DYJetsToLL_M10to50 ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/EWKWMinus2Jets_WToLNu_final.root files_nominal/EWKWMinus2Jets_WToLNu_final.root EWKWMinus2Jets EWKWMinus 0 
./_postAnalyzer_mutau.exe ../files_initial/EWKWPlus2Jets_WToLNu_final.root files_nominal/EWKWPlus2Jets_WToLNu_final.root EWKWPlus2Jets EWKWPlus 0 
./_postAnalyzer_mutau.exe ../files_initial/EWKZ2Jets_ZToLL_final.root files_nominal/EWKZ2Jets_ZToLL_final.root EWKZ2Jets_ZToLL EWKZ2Jets 0 
./_postAnalyzer_mutau.exe ../files_initial/EWKZ2Jets_ZToNuNu_final.root files_nominal/EWKZ2Jets_ZToNuNu_final.root EWKZ2Jets_ZToNuNu EWKZ2Jets 0 
./_postAnalyzer_mutau.exe ../files_initial/GluGluHToTauTau_final.root files_nominal/GluGluHToTauTau_final.root GluGluHToTauTau GluGluH 0 
./_postAnalyzer_mutau.exe ../files_initial/GluGluHToWWTo2L2Nu_final.root files_nominal/GluGluHToWWTo2L2Nu_final.root GluGluHToWWTo2L2Nu GluGluH 0 
./_postAnalyzer_mutau.exe ../files_initial/GluGluZH_HToWW__final.root files_nominal/GluGluZH_HToWW__final.root GluGluZH_HToWW GluGluZH 0 
./_postAnalyzer_mutau.exe ../files_initial/HWminusJ_HToWW_final.root files_nominal/HWminusJ_HToWW_final.root HWminusJ_HToWW HWminusJ 0 
./_postAnalyzer_mutau.exe ../files_initial/HWplusJ_HToWW_final.root files_nominal/HWplusJ_HToWW_final.root HWplusJ_HToWW HWplusJ 0 
./_postAnalyzer_mutau.exe ../files_initial/WWW_final.root files_nominal/WWW_final.root WWW VVV 0 
./_postAnalyzer_mutau.exe ../files_initial/HZJ_HToWW_final.root files_nominal/HZJ_HToWW_final.root HZJ_HToWW HZJ 0 
./_postAnalyzer_mutau.exe ../files_initial/ST_t-channel_antitop_final.root files_nominal/ST_t-channel_antitop_final.root ST_t-channel_antitop ST_t 0 
./_postAnalyzer_mutau.exe ../files_initial/ST_t-channel_top_final.root files_nominal/ST_t-channel_top_final.root ST_t-channel_top ST_t 0 
./_postAnalyzer_mutau.exe ../files_initial/ST_tW_antitop_final.root files_nominal/ST_tW_antitop_final.root ST_tW_antitop ST_t 0 
./_postAnalyzer_mutau.exe ../files_initial/ST_tW_top_final.root files_nominal/ST_tW_top_final.root ST_tW_top ST_t 0 
./_postAnalyzer_mutau.exe ../files_initial/SingleMuonA_final.root files_nominal/SingleMuonA_final.root data_obs data_obs 0 
./_postAnalyzer_mutau.exe ../files_initial/SingleMuonB_final.root files_nominal/SingleMuonB_final.root data_obs data_obs 0 
./_postAnalyzer_mutau.exe ../files_initial/SingleMuonC_final.root files_nominal/SingleMuonC_final.root data_obs data_obs 0 
./_postAnalyzer_mutau.exe ../files_initial/WWZ_final.root files_nominal/WWZ_final.root WWZ VVV 0 
./_postAnalyzer_mutau.exe ../files_initial/WZZ_final.root files_nominal/WZZ_final.root WZZ VVV 0 
./_postAnalyzer_mutau.exe ../files_initial/TTTo2L2Nu_powheg_final.root files_nominal/TTTo2L2Nu_powheg_final.root TTTo2L2Nu TT 0 
./_postAnalyzer_mutau.exe ../files_initial/TTToHadronic_powheg_final.root files_nominal/TTToHadronic_powheg_final.root TTToHadronic TT 0 
./_postAnalyzer_mutau.exe ../files_initial/TTToSemiLeptonic_powheg_final.root files_nominal/TTToSemiLeptonic_powheg_final.root TTToSemiLeptonic TT 0 
./_postAnalyzer_mutau.exe ../files_initial/VBFHToTauTau_final.root files_nominal/VBFHToTauTau_final.root VBFHToTauTau VBFH 0 
./_postAnalyzer_mutau.exe ../files_initial/VBFHToWWTo2L2Nu_final.root files_nominal/VBFHToWWTo2L2Nu_final.root VBFHToWWTo2L2Nu VBFH 0 
./_postAnalyzer_mutau.exe ../files_initial/VVTo2L2Nu_final.root files_nominal/VVTo2L2Nu_final.root VVTo2L2Nu VV 0 
./_postAnalyzer_mutau.exe ../files_initial/W1JetsToLNu_final.root files_nominal/W1JetsToLNu_final.root W1JetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/W2JetsToLNu_final.root files_nominal/W2JetsToLNu_final.root W2JetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/W3JetsToLNu_final.root files_nominal/W3JetsToLNu_final.root W3JetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/W4JetsToLNu_final.root files_nominal/W4JetsToLNu_final.root W4JetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/WGToLNuG_final.root files_nominal/WGToLNuG_final.root WGToLNuG WGToLNuG 0 
./_postAnalyzer_mutau.exe ../files_initial/WJetsToLNu_final.root files_nominal/WJetsToLNu_final.root WJetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/WWTo1L1Nu2Q_final.root files_nominal/WWTo1L1Nu2Q_final.root WWTo1L1Nu2Q VV 0 
./_postAnalyzer_mutau.exe ../files_initial/WWToLNuQQ_final.root files_nominal/WWToLNuQQ_final.root WWToLNuQQ VV 0 
./_postAnalyzer_mutau.exe ../files_initial/WZTo3LNu_final.root files_nominal/WZTo3LNu_final.root WZTo3LNu VV 0 
./_postAnalyzer_mutau.exe ../files_initial/WminusHToTauTau_final.root files_nominal/WminusHToTauTau_final.root WminusHToTauTau WminusH 0 
./_postAnalyzer_mutau.exe ../files_initial/WplusHToTauTau_final.root files_nominal/WplusHToTauTau_final.root WplusHToTauTau WplusH 0 
./_postAnalyzer_mutau.exe ../files_initial/ZHToTauTau_final.root files_nominal/ZHToTauTau_final.root ZHToTauTau ZH 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT100-200_final.root files_nominal/ZJetsToNuNu_HT100-200_final.root ZJetsToNuNu_HT100-200 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT1200-2500_final.root files_nominal/ZJetsToNuNu_HT1200-2500_final.root ZJetsToNuNu_HT1200-2500 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT200-400_final.root files_nominal/ZJetsToNuNu_HT200-400_final.root ZJetsToNuNu_HT200-400 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT2500-Inf_final.root files_nominal/ZJetsToNuNu_HT2500-Inf_final.root ZJetsToNuNu_HT2500-Inf ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT400-600_final.root files_nominal/ZJetsToNuNu_HT400-600_final.root ZJetsToNuNu_HT400-600 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT600-800_final.root files_nominal/ZJetsToNuNu_HT600-800_final.root ZJetsToNuNu_HT600-800 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT800-1200_final.root files_nominal/ZJetsToNuNu_HT800-1200_final.root ZJetsToNuNu_HT800-1200 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZZTo2L2Q_final.root files_nominal/ZZTo2L2Q_final.root ZZTo2L2Q VV 0 
./_postAnalyzer_mutau.exe ../files_initial/ZZTo4L_final.root files_nominal/ZZTo4L_final.root ZZTo4L VV 0 
./_postAnalyzer_mutau.exe ../files_initial/ZZZ_final.root files_nominal/ZZZ_final.root ZZZ VVV 0 
./_postAnalyzer_mutau.exe ../files_initial/ggZH_HToTauTau_ZToLL_final.root files_nominal/ggZH_HToTauTau_ZToLL_final.root ggZH_HToTauTau_ZToLL ggZH 0 
./_postAnalyzer_mutau.exe ../files_initial/ggZH_HToTauTau_ZToNuNu_final.root files_nominal/ggZH_HToTauTau_ZToNuNu_final.root ggZH_HToTauTau_ZToNuNu ggZH 0 
./_postAnalyzer_mutau.exe ../files_initial/ggZH_HToTauTau_ZToQQ_final.root files_nominal/ggZH_HToTauTau_ZToQQ_final.root ggZH_HToTauTau_ZToQQ ggZH 0 
./_postAnalyzer_mutau.exe ../files_initial/ttHToNonbb_final.root files_nominal/ttHToNonbb_final.root ttHToNonbb ttH 0 
./_postAnalyzer_mutau.exe ../files_initial/SingleMuonD_PromptReco_final.root files_nominal/SingleMuonD_PromptReco_final.root data_obs data_obs 0 

hadd -f f_mutau_initial.root files_nominal/*.root 
echo "*************** root file made ***************" 
#sh do_plots_mutau.sh 
echo "*************** plots made ***************" 
