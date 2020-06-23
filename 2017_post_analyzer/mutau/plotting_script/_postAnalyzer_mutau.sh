
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

./_postAnalyzer_mutau.exe ../files_initial/DY1JetsToLL_M-50_TuneCP5_final.root files_nominal/DY1JetsToLL_M-50_TuneCP5_final.root DY1JetsToLL ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DY1JetsToLL_M-50_TuneCP5_final.root files_nominal/DY1JetsToLL_M-50_TuneCP5_stitch_final.root ZTT1jet ZTTjet 0 
./_postAnalyzer_mutau.exe ../files_initial/DY2JetsToLL_M-50_TuneCP5_final.root files_nominal/DY2JetsToLL_M-50_TuneCP5_final.root DY2JetsToLL ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DY2JetsToLL_M-50_TuneCP5_final.root files_nominal/DY2JetsToLL_M-50_TuneCP5_stitch_final.root ZTT2jet ZTTjet 0 

./_postAnalyzer_mutau.exe ../files_initial/DY3JetsToLL_M-50_TuneCP5_ext1_final.root files_nominal/DY3JetsToLL_M-50_TuneCP5_ext1_final.root DY3JetsToLL ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DY3JetsToLL_M-50_TuneCP5_ext1_final.root files_nominal/DY3JetsToLL_M-50_TuneCP5_ext1_stitch_final.root ZTT3jet ZTTjet 0 

./_postAnalyzer_mutau.exe ../files_initial/DY3JetsToLL_M-50_TuneCP5_v1_final.root files_nominal/DY3JetsToLL_M-50_TuneCP5_v1_final.root DY3JetsToLL ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DY3JetsToLL_M-50_TuneCP5_v1_final.root files_nominal/DY3JetsToLL_M-50_TuneCP5_v1_stitch_final.root ZTT3jet ZTTjet 0 

./_postAnalyzer_mutau.exe ../files_initial/DY4JetsToLL_M-50_TuneCP5_final.root files_nominal/DY4JetsToLL_M-50_TuneCP5_final.root DY4JetsToLL ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DY4JetsToLL_M-50_TuneCP5_final.root files_nominal/DY4JetsToLL_M-50_TuneCP5_stitch_final.root ZTT4jet ZTTjet 0 
#./_postAnalyzer_mutau.exe ../files_initial/DYJetsToLL_M-10to50_TuneCP5_final.root files_nominal/DYJetsToLL_M-10to50_TuneCP5_final.root DYJetsToLL ZTT 0 

./_postAnalyzer_mutau.exe ../files_initial/DYJetsToLL_M-50_TuneCP5_ext1_v1_final.root files_nominal/DYJetsToLL_M-50_TuneCP5_ext1_v1_final.root DYJetsToLL ZTT 0 
./_postAnalyzer_mutau.exe ../files_initial/DYJetsToLL_M-50_TuneCP5_ext1_v1_final.root files_nominal/DYJetsToLL_M-50_TuneCP5_ext1_v1_stitch_final.root ZTTjet_inc ZTTjet 0 

./_postAnalyzer_mutau.exe ../files_initial/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_final.root files_nominal/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_final.root EWKWMinus2Jets EWKWMinus 0 
./_postAnalyzer_mutau.exe ../files_initial/EWKWPlus2Jets_WToLNu_M-50_TuneCP5_final.root files_nominal/EWKWPlus2Jets_WToLNu_M-50_TuneCP5_final.root EWKWPlus2Jets EWKWPlus 0 
./_postAnalyzer_mutau.exe ../files_initial/EWKZ2Jets_ZToLL_M-50_TuneCP5_final.root files_nominal/EWKZ2Jets_ZToLL_M-50_TuneCP5_final.root EWKZ2Jets_ZToLL EWKZ2Jets 0 
./_postAnalyzer_mutau.exe ../files_initial/EWKZ2Jets_ZToNuNu_TuneCP5_final.root files_nominal/EWKZ2Jets_ZToNuNu_TuneCP5_final.root EWKZ2Jets_ZToNuNu EWKZ2Jets 0 
./_postAnalyzer_mutau.exe ../files_initial/GluGluHToTauTau_M125_final.root files_nominal/GluGluHToTauTau_M125_final.root GluGluHToTauTau GluGluH 0 
./_postAnalyzer_mutau.exe ../files_initial/GluGluHToWWTo2L2Nu_M125_final.root files_nominal/GluGluHToWWTo2L2Nu_M125_final.root GluGluHToWWTo2L2Nu GluGluH 0 
./_postAnalyzer_mutau.exe ../files_initial/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_final.root files_nominal/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_final.root ST_t-channel_antitop ST_t 0 
./_postAnalyzer_mutau.exe ../files_initial/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_final.root files_nominal/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_final.root ST_t-channel_top ST_t 0 
./_postAnalyzer_mutau.exe ../files_initial/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_final.root files_nominal/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_final.root ST_tW_antitop ST_t 0 
./_postAnalyzer_mutau.exe ../files_initial/ST_tW_top_5f_inclusiveDecays_TuneCP5_final.root files_nominal/ST_tW_top_5f_inclusiveDecays_TuneCP5_final.root ST_tW_top ST_t 0 
./_postAnalyzer_mutau.exe ../files_initial/SingleMuon_EraB_final.root files_nominal/SingleMuon_EraB_final.root data_obs data_obs 0 
./_postAnalyzer_mutau.exe ../files_initial/SingleMuon_EraC_final.root files_nominal/SingleMuon_EraC_final.root data_obs data_obs 0 
./_postAnalyzer_mutau.exe ../files_initial/SingleMuon_EraD_final.root files_nominal/SingleMuon_EraD_final.root data_obs data_obs 0 
./_postAnalyzer_mutau.exe ../files_initial/SingleMuon_EraE_final.root files_nominal/SingleMuon_EraE_final.root data_obs data_obs 0 
./_postAnalyzer_mutau.exe ../files_initial/SingleMuon_EraF_final.root files_nominal/SingleMuon_EraF_final.root data_obs data_obs 0 
./_postAnalyzer_mutau.exe ../files_initial/TTTo2L2Nu_TuneCP5_final.root files_nominal/TTTo2L2Nu_TuneCP5_final.root TTTo2L2Nu TT 0 
./_postAnalyzer_mutau.exe ../files_initial/TTToHadronic_TuneCP5_final.root files_nominal/TTToHadronic_TuneCP5_final.root TTToHadronic TT 0 
./_postAnalyzer_mutau.exe ../files_initial/TTToSemiLeptonic_TuneCP5_final.root files_nominal/TTToSemiLeptonic_TuneCP5_final.root TTToSemiLeptonic TT 0 
./_postAnalyzer_mutau.exe ../files_initial/VBFHToTauTau_M125_final.root files_nominal/VBFHToTauTau_M125_final.root VBFHToTauTau VBFH 0 
./_postAnalyzer_mutau.exe ../files_initial/VBFHToWWTo2L2Nu_M125_final.root files_nominal/VBFHToWWTo2L2Nu_M125_final.root VBFHToWWTo2L2Nu VBFH 0 
./_postAnalyzer_mutau.exe ../files_initial/VVTo2L2Nu_final.root files_nominal/VVTo2L2Nu_final.root VVTo2L2Nu VV 0 
./_postAnalyzer_mutau.exe ../files_initial/W1JetsToLNu_TuneCP5_final.root files_nominal/W1JetsToLNu_TuneCP5_final.root W1JetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/W1JetsToLNu_TuneCP5_final.root files_nominal/W1JetsToLNu_TuneCP5_stitch_final.root W1Jet WJets_jets 0 
./_postAnalyzer_mutau.exe ../files_initial/W2JetsToLNu_TuneCP5_final.root files_nominal/W2JetsToLNu_TuneCP5_final.root W2JetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/W2JetsToLNu_TuneCP5_final.root files_nominal/W2JetsToLNu_TuneCP5_stitch_final.root W2Jet WJets_jets 0 
./_postAnalyzer_mutau.exe ../files_initial/W3JetsToLNu_TuneCP5_final.root files_nominal/W3JetsToLNu_TuneCP5_final.root W3JetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/W3JetsToLNu_TuneCP5_final.root files_nominal/W3JetsToLNu_TuneCP5_stitch_final.root W3Jet WJets_jets 0 
./_postAnalyzer_mutau.exe ../files_initial/W4JetsToLNu_TuneCP5_final.root files_nominal/W4JetsToLNu_TuneCP5_final.root W4JetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/W4JetsToLNu_TuneCP5_final.root files_nominal/W4JetsToLNu_TuneCP5_stitch_final.root W4Jet WJets_jets 0 
./_postAnalyzer_mutau.exe ../files_initial/WJetsToLNu_TuneCP5_final.root files_nominal/WJetsToLNu_TuneCP5_final.root WJetsToLNu WJets 0 
./_postAnalyzer_mutau.exe ../files_initial/WJetsToLNu_TuneCP5_final.root files_nominal/WJetsToLNu_TuneCP5_stitch_final.root WJets_inc WJets_jets 0 
./_postAnalyzer_mutau.exe ../files_initial/WWTo1L1Nu2Q_final.root files_nominal/WWTo1L1Nu2Q_final.root WWTo1L1Nu2Q VV 0 
./_postAnalyzer_mutau.exe ../files_initial/WWToLNuQQ_NNPDF31_TuneCP5_final.root files_nominal/WWToLNuQQ_NNPDF31_TuneCP5_final.root WWToLNuQQ VV 0 
./_postAnalyzer_mutau.exe ../files_initial/WWW_4F_TuneCP5_final.root files_nominal/WWW_4F_TuneCP5_final.root WWW VVV 0 
./_postAnalyzer_mutau.exe ../files_initial/WWZ_4F_TuneCP5_final.root files_nominal/WWZ_4F_TuneCP5_final.root WWZ VVV 0 
./_postAnalyzer_mutau.exe ../files_initial/WW_TuneCP5_final.root files_nominal/WW_TuneCP5_final.root WW VV 0 
./_postAnalyzer_mutau.exe ../files_initial/WZTo3LNu_TuneCP5_final.root files_nominal/WZTo3LNu_TuneCP5_final.root WZTo3LNu VV 0 
./_postAnalyzer_mutau.exe ../files_initial/WZZ_TuneCP5_final.root files_nominal/WZZ_TuneCP5_final.root WZZ VVV 0 
./_postAnalyzer_mutau.exe ../files_initial/WZ_TuneCP5_final.root files_nominal/WZ_TuneCP5_final.root WZ VV 0 
./_postAnalyzer_mutau.exe ../files_initial/WminusHToTauTau_M125_final.root files_nominal/WminusHToTauTau_M125_final.root WminusHToTauTau WminusH 0 
./_postAnalyzer_mutau.exe ../files_initial/WplusHToTauTau_M125_final.root files_nominal/WplusHToTauTau_M125_final.root WplusHToTauTau WplusH 0 
./_postAnalyzer_mutau.exe ../files_initial/ZHToTauTau_M125_final.root files_nominal/ZHToTauTau_M125_final.root ZHToTauTau ZH 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT-100To200_final.root files_nominal/ZJetsToNuNu_HT-100To200_final.root ZJetsToNuNu_HT100-200 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT-1200To2500_final.root files_nominal/ZJetsToNuNu_HT-1200To2500_final.root ZJetsToNuNu_HT1200-2500 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT-200To400_final.root files_nominal/ZJetsToNuNu_HT-200To400_final.root ZJetsToNuNu_HT200-400 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT-2500ToInf_final.root files_nominal/ZJetsToNuNu_HT-2500ToInf_final.root ZJetsToNuNu_HT2500-Inf ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT-400To600_final.root files_nominal/ZJetsToNuNu_HT-400To600_final.root ZJetsToNuNu_HT400-600 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT-600To800_final.root files_nominal/ZJetsToNuNu_HT-600To800_final.root ZJetsToNuNu_HT600-800 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZJetsToNuNu_HT-800To1200_final.root files_nominal/ZJetsToNuNu_HT-800To1200_final.root ZJetsToNuNu_HT800-1200 ZJetsToNuNu 0 
./_postAnalyzer_mutau.exe ../files_initial/ZZTo2L2Q_final.root files_nominal/ZZTo2L2Q_final.root ZZTo2L2Q VV 0 
./_postAnalyzer_mutau.exe ../files_initial/ZZTo4L_TuneCP5_final.root files_nominal/ZZTo4L_TuneCP5_final.root ZZTo4L VV 0 
./_postAnalyzer_mutau.exe ../files_initial/ZZZ_TuneCP5_final.root files_nominal/ZZZ_TuneCP5_final.root ZZZ VVV 0 
./_postAnalyzer_mutau.exe ../files_initial/ZZ_TuneCP5_final.root files_nominal/ZZ_TuneCP5_final.root ZZ VV 0 

hadd -f f_mutau_initial.root files_nominal/*.root 
echo "*************** root file made ***************" 
sh do_make_plots.sh
echo "*************** plots made ***************" 
