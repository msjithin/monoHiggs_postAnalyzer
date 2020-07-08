
#!/bin/bash
set -e 
outDIR="files_initial"
if [ -d "$outDIR" ]; then
 echo "$outDIR exists"
 if [ "$(ls -A $outDIR)" ]; then
 echo "Take action $outDIR is not Empty .... removing existing files ....."
 rm $outDIR/*.root
 else
 echo " $outDIR is Empty"
 fi
else
 echo "$outDIR created"
 mkdir $outDIR
fi
hadd $outDIR/DY1JetsToLL_M-50_TuneCP5_final.root rootFiles/DY1JetsToLL_M-50_TuneCP5_00.root rootFiles/DY1JetsToLL_M-50_TuneCP5_01.root rootFiles/DY1JetsToLL_M-50_TuneCP5_02.root rootFiles/DY1JetsToLL_M-50_TuneCP5_03.root rootFiles/DY1JetsToLL_M-50_TuneCP5_04.root 
hadd $outDIR/DY2JetsToLL_M-50_TuneCP5_final.root rootFiles/DY2JetsToLL_M-50_TuneCP5_00.root rootFiles/DY2JetsToLL_M-50_TuneCP5_01.root 
hadd $outDIR/DY3JetsToLL_M-50_TuneCP5_ext1_final.root rootFiles/DY3JetsToLL_M-50_TuneCP5_ext1_00.root rootFiles/DY3JetsToLL_M-50_TuneCP5_v1_00.root 
hadd $outDIR/DY4JetsToLL_M-50_TuneCP5_final.root rootFiles/DY4JetsToLL_M-50_TuneCP5_00.root 
hadd $outDIR/DYJetsToLL_M-10to50_TuneCP5_final.root rootFiles/DYJetsToLL_M-10to50_TuneCP5_00.root rootFiles/DYJetsToLL_M-10to50_TuneCP5_01.root rootFiles/DYJetsToLL_M-10to50_TuneCP5_02.root rootFiles/DYJetsToLL_M-10to50_TuneCP5_03.root rootFiles/DYJetsToLL_M-10to50_TuneCP5_04.root 
hadd $outDIR/DYJetsToLL_M-50_TuneCP5_ext1_v1_final.root rootFiles/DYJetsToLL_M-50_TuneCP5_ext1_v1_00.root rootFiles/DYJetsToLL_M-50_TuneCP5_ext1_v1_01.root rootFiles/DYJetsToLL_M-50_TuneCP5_ext1_v1_02.root rootFiles/DYJetsToLL_M-50_TuneCP5_ext1_v1_03.root rootFiles/DYJetsToLL_M-50_TuneCP5_ext1_v1_04.root rootFiles/DYJetsToLL_M-50_TuneCP5_ext1_v1_05.root rootFiles/DYJetsToLL_M-50_TuneCP5_v1_00.root rootFiles/DYJetsToLL_M-50_TuneCP5_v1_01.root rootFiles/DYJetsToLL_M-50_TuneCP5_v1_02.root rootFiles/DYJetsToLL_M-50_TuneCP5_v1_03.root rootFiles/DYJetsToLL_M-50_TuneCP5_v1_04.root rootFiles/DYJetsToLL_M-50_TuneCP5_v1_05.root 
hadd $outDIR/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_final.root rootFiles/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_00.root 
hadd $outDIR/EWKWPlus2Jets_WToLNu_M-50_TuneCP5_final.root rootFiles/EWKWPlus2Jets_WToLNu_M-50_TuneCP5_00.root 
hadd $outDIR/EWKZ2Jets_ZToLL_M-50_TuneCP5_final.root rootFiles/EWKZ2Jets_ZToLL_M-50_TuneCP5_00.root 
hadd $outDIR/EWKZ2Jets_ZToNuNu_TuneCP5_final.root rootFiles/EWKZ2Jets_ZToNuNu_TuneCP5_00.root 
hadd $outDIR/GluGluHToTauTau_M125_final.root rootFiles/GluGluHToTauTau_M125_00.root 
hadd $outDIR/GluGluHToWWTo2L2Nu_M125_final.root rootFiles/GluGluHToWWTo2L2Nu_M125_00.root 
hadd $outDIR/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_final.root rootFiles/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_00.root 
hadd $outDIR/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_final.root rootFiles/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_00.root 
hadd $outDIR/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_final.root rootFiles/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_00.root 
hadd $outDIR/ST_tW_top_5f_inclusiveDecays_TuneCP5_final.root rootFiles/ST_tW_top_5f_inclusiveDecays_TuneCP5_00.root 
hadd $outDIR/TTTo2L2Nu_TuneCP5_final.root rootFiles/TTTo2L2Nu_TuneCP5_00.root 
hadd $outDIR/TTToHadronic_TuneCP5_final.root rootFiles/TTToHadronic_TuneCP5_00.root rootFiles/TTToHadronic_TuneCP5_01.root rootFiles/TTToHadronic_TuneCP5_02.root rootFiles/TTToHadronic_TuneCP5_03.root rootFiles/TTToHadronic_TuneCP5_04.root rootFiles/TTToHadronic_TuneCP5_05.root rootFiles/TTToHadronic_TuneCP5_06.root rootFiles/TTToHadronic_TuneCP5_07.root rootFiles/TTToHadronic_TuneCP5_08.root 
hadd $outDIR/TTToSemiLeptonic_TuneCP5_final.root rootFiles/TTToSemiLeptonic_TuneCP5_00.root rootFiles/TTToSemiLeptonic_TuneCP5_01.root rootFiles/TTToSemiLeptonic_TuneCP5_02.root rootFiles/TTToSemiLeptonic_TuneCP5_03.root rootFiles/TTToSemiLeptonic_TuneCP5_04.root 
hadd $outDIR/Tau_EraB_final.root rootFiles/Tau_EraB_00.root rootFiles/Tau_EraB_01.root 
hadd $outDIR/Tau_EraC_final.root rootFiles/Tau_EraC_00.root rootFiles/Tau_EraC_01.root rootFiles/Tau_EraC_02.root rootFiles/Tau_EraC_03.root 
hadd $outDIR/Tau_EraD_final.root rootFiles/Tau_EraD_00.root rootFiles/Tau_EraD_01.root rootFiles/Tau_EraD_02.root 
hadd $outDIR/Tau_EraE_final.root rootFiles/Tau_EraE_00.root rootFiles/Tau_EraE_01.root rootFiles/Tau_EraE_02.root rootFiles/Tau_EraE_03.root 
hadd $outDIR/Tau_EraF_final.root rootFiles/Tau_EraF_00.root rootFiles/Tau_EraF_01.root rootFiles/Tau_EraF_02.root rootFiles/Tau_EraF_03.root rootFiles/Tau_EraF_04.root 
hadd $outDIR/VBFHToTauTau_M125_final.root rootFiles/VBFHToTauTau_M125_00.root 
hadd $outDIR/VBFHToWWTo2L2Nu_M125_final.root rootFiles/VBFHToWWTo2L2Nu_M125_00.root 
hadd $outDIR/VVTo2L2Nu_final.root rootFiles/VVTo2L2Nu_00.root 
hadd $outDIR/W1JetsToLNu_TuneCP5_final.root rootFiles/W1JetsToLNu_TuneCP5_00.root rootFiles/W1JetsToLNu_TuneCP5_01.root rootFiles/W1JetsToLNu_TuneCP5_02.root rootFiles/W1JetsToLNu_TuneCP5_03.root rootFiles/W1JetsToLNu_TuneCP5_04.root rootFiles/W1JetsToLNu_TuneCP5_05.root rootFiles/W1JetsToLNu_TuneCP5_06.root 
hadd $outDIR/W2JetsToLNu_TuneCP5_final.root rootFiles/W2JetsToLNu_TuneCP5_00.root 
hadd $outDIR/W3JetsToLNu_TuneCP5_final.root rootFiles/W3JetsToLNu_TuneCP5_00.root rootFiles/W3JetsToLNu_TuneCP5_01.root 
hadd $outDIR/W4JetsToLNu_TuneCP5_final.root rootFiles/W4JetsToLNu_TuneCP5_00.root rootFiles/W4JetsToLNu_TuneCP5_01.root 
hadd $outDIR/WJetsToLNu_TuneCP5_final.root rootFiles/WJetsToLNu_TuneCP5_00.root rootFiles/WJetsToLNu_TuneCP5_01.root rootFiles/WJetsToLNu_TuneCP5_02.root rootFiles/WJetsToLNu_TuneCP5_03.root rootFiles/WJetsToLNu_TuneCP5_04.root 
hadd $outDIR/WWTo1L1Nu2Q_final.root rootFiles/WWTo1L1Nu2Q_00.root 
hadd $outDIR/WWToLNuQQ_NNPDF31_TuneCP5_final.root rootFiles/WWToLNuQQ_NNPDF31_TuneCP5_00.root 
hadd $outDIR/WWW_4F_TuneCP5_final.root rootFiles/WWW_4F_TuneCP5_00.root 
hadd $outDIR/WWZ_4F_TuneCP5_final.root rootFiles/WWZ_4F_TuneCP5_00.root 
hadd $outDIR/WW_TuneCP5_final.root rootFiles/WW_TuneCP5_00.root 
hadd $outDIR/WZTo3LNu_TuneCP5_final.root rootFiles/WZTo3LNu_TuneCP5_00.root rootFiles/WZTo3LNu_TuneCP5_01.root 
hadd $outDIR/WZZ_TuneCP5_final.root rootFiles/WZZ_TuneCP5_00.root 
hadd $outDIR/WZ_TuneCP5_final.root rootFiles/WZ_TuneCP5_00.root 
hadd $outDIR/WminusHToTauTau_M125_final.root rootFiles/WminusHToTauTau_M125_00.root 
hadd $outDIR/WplusHToTauTau_M125_final.root rootFiles/WplusHToTauTau_M125_00.root 
hadd $outDIR/ZHToTauTau_M125_final.root rootFiles/ZHToTauTau_M125_00.root 
hadd $outDIR/ZJetsToNuNu_HT-100To200_final.root rootFiles/ZJetsToNuNu_HT-100To200_00.root rootFiles/ZJetsToNuNu_HT-100To200_01.root rootFiles/ZJetsToNuNu_HT-100To200_02.root 
hadd $outDIR/ZJetsToNuNu_HT-1200To2500_final.root rootFiles/ZJetsToNuNu_HT-1200To2500_00.root 
hadd $outDIR/ZJetsToNuNu_HT-200To400_final.root rootFiles/ZJetsToNuNu_HT-200To400_00.root rootFiles/ZJetsToNuNu_HT-200To400_01.root rootFiles/ZJetsToNuNu_HT-200To400_02.root 
hadd $outDIR/ZJetsToNuNu_HT-2500ToInf_final.root rootFiles/ZJetsToNuNu_HT-2500ToInf_00.root 
hadd $outDIR/ZJetsToNuNu_HT-400To600_final.root rootFiles/ZJetsToNuNu_HT-400To600_00.root 
hadd $outDIR/ZJetsToNuNu_HT-600To800_final.root rootFiles/ZJetsToNuNu_HT-600To800_00.root 
hadd $outDIR/ZJetsToNuNu_HT-800To1200_final.root rootFiles/ZJetsToNuNu_HT-800To1200_00.root 
hadd $outDIR/ZZTo2L2Q_final.root rootFiles/ZZTo2L2Q_00.root rootFiles/ZZTo2L2Q_01.root rootFiles/ZZTo2L2Q_02.root 
hadd $outDIR/ZZTo4L_TuneCP5_final.root rootFiles/ZZTo4L_TuneCP5_00.root rootFiles/ZZTo4L_TuneCP5_01.root rootFiles/ZZTo4L_TuneCP5_02.root 
hadd $outDIR/ZZZ_TuneCP5_final.root rootFiles/ZZZ_TuneCP5_00.root 
hadd $outDIR/ZZ_TuneCP5_final.root rootFiles/ZZ_TuneCP5_00.root 


