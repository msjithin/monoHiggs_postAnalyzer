
#!/bin/bash
set -e 

if [ "$(ls -A files_initial)" ]; then
echo "Take action files_initial/ is not Empty .... removing existing files ....."
rm files_initial/*.root
else
echo " files_initial/ is Empty"
fi

hadd files_initial/DY1JetsToLL_final.root rootFiles/DY1JetsToLL_00.root rootFiles/DY1JetsToLL_01.root rootFiles/DY1JetsToLL_02.root rootFiles/DY1JetsToLL_03.root rootFiles/DY1JetsToLL_04.root 
hadd files_initial/DY2JetsToLL_final.root rootFiles/DY2JetsToLL_00.root rootFiles/DY2JetsToLL_01.root rootFiles/DY2JetsToLL_02.root 
hadd files_initial/DY3JetsToLL_final.root rootFiles/DY3JetsToLL_00.root 
hadd files_initial/DY4JetsToLL_final.root rootFiles/DY4JetsToLL_00.root 
hadd files_initial/DYJetsToLL_final.root rootFiles/DYJetsToLL_00.root rootFiles/DYJetsToLL_01.root rootFiles/DYJetsToLL_02.root rootFiles/DYJetsToLL_03.root rootFiles/DYJetsToLL_04.root rootFiles/DYJetsToLL_05.root rootFiles/DYJetsToLL_06.root 
hadd files_initial/DYJetsToLL_Incl_HT_final.root rootFiles/DYJetsToLL_Incl_HT_00.root rootFiles/DYJetsToLL_Incl_HT_01.root 
hadd files_initial/DYJetsToLL_M10to50_final.root rootFiles/DYJetsToLL_M10to50_00.root rootFiles/DYJetsToLL_M10to50_01.root rootFiles/DYJetsToLL_M10to50_02.root rootFiles/DYJetsToLL_M10to50_03.root 
hadd files_initial/EGamma2018A_final.root rootFiles/EGamma2018A_00.root rootFiles/EGamma2018A_01.root rootFiles/EGamma2018A_02.root rootFiles/EGamma2018A_03.root rootFiles/EGamma2018A_04.root rootFiles/EGamma2018A_05.root rootFiles/EGamma2018A_06.root rootFiles/EGamma2018A_07.root rootFiles/EGamma2018A_08.root rootFiles/EGamma2018A_09.root rootFiles/EGamma2018A_10.root 
hadd files_initial/WWW_final.root rootFiles/WWW_00.root 
hadd files_initial/EGamma2018B_final.root rootFiles/EGamma2018B_00.root rootFiles/EGamma2018B_01.root rootFiles/EGamma2018B_02.root rootFiles/EGamma2018B_03.root rootFiles/EGamma2018B_04.root rootFiles/EGamma2018B_05.root 
hadd files_initial/EGamma2018C_final.root rootFiles/EGamma2018C_00.root rootFiles/EGamma2018C_01.root rootFiles/EGamma2018C_02.root rootFiles/EGamma2018C_03.root rootFiles/EGamma2018C_04.root rootFiles/EGamma2018C_05.root rootFiles/EGamma2018C_06.root rootFiles/EGamma2018C_07.root 
hadd files_initial/EGamma2018D_final.root rootFiles/EGamma2018D_00.root rootFiles/EGamma2018D_01.root rootFiles/EGamma2018D_02.root rootFiles/EGamma2018D_03.root rootFiles/EGamma2018D_04.root rootFiles/EGamma2018D_05.root rootFiles/EGamma2018D_06.root rootFiles/EGamma2018D_07.root rootFiles/EGamma2018D_08.root rootFiles/EGamma2018D_09.root rootFiles/EGamma2018D_10.root rootFiles/EGamma2018D_11.root rootFiles/EGamma2018D_12.root rootFiles/EGamma2018D_13.root rootFiles/EGamma2018D_14.root rootFiles/EGamma2018D_15.root rootFiles/EGamma2018D_16.root rootFiles/EGamma2018D_17.root rootFiles/EGamma2018D_18.root rootFiles/EGamma2018D_19.root 
hadd files_initial/WWZ_final.root rootFiles/WWZ_00.root 
hadd files_initial/EWKWMinus2Jets_WToLNu_final.root rootFiles/EWKWMinus2Jets_WToLNu_00.root 
hadd files_initial/EWKWPlus2Jets_WToLNu_final.root rootFiles/EWKWPlus2Jets_WToLNu_00.root 
hadd files_initial/EWKZ2Jets_ZToNuNu_final.root rootFiles/EWKZ2Jets_ZToNuNu_00.root 
hadd files_initial/GluGluHToTauTau_final.root rootFiles/GluGluHToTauTau_00.root rootFiles/GluGluHToTauTau_01.root 
hadd files_initial/GluGluZH_HToWW_final.root rootFiles/GluGluZH_HToWW_00.root 
hadd files_initial/HWminusJ_HToWW_final.root rootFiles/HWminusJ_HToWW_00.root 
hadd files_initial/HWplusJ_HToWW_final.root rootFiles/HWplusJ_HToWW_00.root 
hadd files_initial/HZJ_HToWW_final.root rootFiles/HZJ_HToWW_00.root 
hadd files_initial/ST_t-channel_antitop_final.root rootFiles/ST_t-channel_antitop_00.root rootFiles/ST_t-channel_antitop_01.root rootFiles/ST_t-channel_antitop_02.root rootFiles/ST_t-channel_antitop_03.root rootFiles/ST_t-channel_antitop_04.root rootFiles/ST_t-channel_antitop_05.root rootFiles/ST_t-channel_antitop_06.root rootFiles/ST_t-channel_antitop_07.root 
hadd files_initial/ST_t-channel_top_final.root rootFiles/ST_t-channel_top_00.root rootFiles/ST_t-channel_top_01.root rootFiles/ST_t-channel_top_02.root rootFiles/ST_t-channel_top_03.root rootFiles/ST_t-channel_top_04.root rootFiles/ST_t-channel_top_05.root rootFiles/ST_t-channel_top_06.root rootFiles/ST_t-channel_top_07.root rootFiles/ST_t-channel_top_08.root rootFiles/ST_t-channel_top_09.root rootFiles/ST_t-channel_top_10.root rootFiles/ST_t-channel_top_11.root 
hadd files_initial/WZZ_final.root rootFiles/WZZ_00.root 
hadd files_initial/ST_tW_antitop_final.root rootFiles/ST_tW_antitop_00.root 
hadd files_initial/ST_tW_top_final.root rootFiles/ST_tW_top_00.root 
hadd files_initial/TTTo2L2Nu_powheg_final.root rootFiles/TTTo2L2Nu_powheg_00.root rootFiles/TTTo2L2Nu_powheg_01.root rootFiles/TTTo2L2Nu_powheg_02.root rootFiles/TTTo2L2Nu_powheg_03.root rootFiles/TTTo2L2Nu_powheg_04.root rootFiles/TTTo2L2Nu_powheg_05.root rootFiles/TTTo2L2Nu_powheg_06.root 
hadd files_initial/TTToHadronic_powheg_final.root rootFiles/TTToHadronic_powheg_00.root rootFiles/TTToHadronic_powheg_01.root rootFiles/TTToHadronic_powheg_02.root rootFiles/TTToHadronic_powheg_03.root rootFiles/TTToHadronic_powheg_04.root rootFiles/TTToHadronic_powheg_05.root rootFiles/TTToHadronic_powheg_06.root rootFiles/TTToHadronic_powheg_07.root rootFiles/TTToHadronic_powheg_08.root rootFiles/TTToHadronic_powheg_09.root rootFiles/TTToHadronic_powheg_10.root rootFiles/TTToHadronic_powheg_11.root rootFiles/TTToHadronic_powheg_12.root rootFiles/TTToHadronic_powheg_13.root rootFiles/TTToHadronic_powheg_14.root 
hadd files_initial/ZZTo4L_final.root rootFiles/ZZTo4L_00.root 
hadd files_initial/TTToSemiLeptonic_powheg_final.root rootFiles/TTToSemiLeptonic_powheg_00.root rootFiles/TTToSemiLeptonic_powheg_01.root rootFiles/TTToSemiLeptonic_powheg_02.root rootFiles/TTToSemiLeptonic_powheg_03.root rootFiles/TTToSemiLeptonic_powheg_04.root rootFiles/TTToSemiLeptonic_powheg_05.root rootFiles/TTToSemiLeptonic_powheg_06.root rootFiles/TTToSemiLeptonic_powheg_07.root rootFiles/TTToSemiLeptonic_powheg_08.root rootFiles/TTToSemiLeptonic_powheg_09.root rootFiles/TTToSemiLeptonic_powheg_10.root rootFiles/TTToSemiLeptonic_powheg_11.root rootFiles/TTToSemiLeptonic_powheg_12.root rootFiles/TTToSemiLeptonic_powheg_13.root rootFiles/TTToSemiLeptonic_powheg_14.root rootFiles/TTToSemiLeptonic_powheg_15.root 
hadd files_initial/VBFHToTauTau_final.root rootFiles/VBFHToTauTau_00.root 
hadd files_initial/VBFHToWWTo2L2Nu_final.root rootFiles/VBFHToWWTo2L2Nu_00.root 
hadd files_initial/W1JetsToLNu_final.root rootFiles/W1JetsToLNu_00.root rootFiles/W1JetsToLNu_01.root rootFiles/W1JetsToLNu_02.root rootFiles/W1JetsToLNu_03.root rootFiles/W1JetsToLNu_04.root rootFiles/W1JetsToLNu_05.root 
hadd files_initial/W2JetsToLNu_final.root rootFiles/W2JetsToLNu_00.root rootFiles/W2JetsToLNu_01.root rootFiles/W2JetsToLNu_02.root 
hadd files_initial/ZZZ_final.root rootFiles/ZZZ_00.root 
hadd files_initial/W3JetsToLNu_final.root rootFiles/W3JetsToLNu_00.root rootFiles/W3JetsToLNu_01.root 
hadd files_initial/W4JetsToLNu_final.root rootFiles/W4JetsToLNu_00.root rootFiles/W4JetsToLNu_01.root 
hadd files_initial/WGToLNuG_final.root rootFiles/WGToLNuG_00.root 
hadd files_initial/WJetsToLNu_Incl_final.root rootFiles/WJetsToLNu_Incl_00.root rootFiles/WJetsToLNu_Incl_01.root rootFiles/WJetsToLNu_Incl_02.root rootFiles/WJetsToLNu_Incl_03.root rootFiles/WJetsToLNu_Incl_04.root rootFiles/WJetsToLNu_Incl_05.root rootFiles/WJetsToLNu_Incl_06.root rootFiles/WJetsToLNu_Incl_07.root rootFiles/WJetsToLNu_Incl_08.root 
hadd files_initial/WJetsToLNu_Incl_HT_final.root rootFiles/WJetsToLNu_Incl_HT_00.root rootFiles/WJetsToLNu_Incl_HT_01.root rootFiles/WJetsToLNu_Incl_HT_02.root rootFiles/WJetsToLNu_Incl_HT_03.root rootFiles/WJetsToLNu_Incl_HT_04.root rootFiles/WJetsToLNu_Incl_HT_05.root rootFiles/WJetsToLNu_Incl_HT_06.root rootFiles/WJetsToLNu_Incl_HT_07.root rootFiles/WJetsToLNu_Incl_HT_08.root 
hadd files_initial/WWTo1L1Nu2Q_final.root rootFiles/WWTo1L1Nu2Q_00.root 
hadd files_initial/WWToLNuQQ_final.root rootFiles/WWToLNuQQ_00.root rootFiles/WWToLNuQQ_01.root 
hadd files_initial/WZTo3LNu_final.root rootFiles/WZTo3LNu_00.root rootFiles/WZTo3LNu_01.root 
hadd files_initial/WplusHToTauTau_final.root rootFiles/WplusHToTauTau_00.root 
hadd files_initial/ZHToTauTau_final.root rootFiles/ZHToTauTau_00.root 
hadd files_initial/ZJetsToNuNu_HT1200-2500_final.root rootFiles/ZJetsToNuNu_HT1200-2500_00.root 
hadd files_initial/ZJetsToNuNu_HT200-400_final.root rootFiles/ZJetsToNuNu_HT200-400_00.root rootFiles/ZJetsToNuNu_HT200-400_01.root rootFiles/ZJetsToNuNu_HT200-400_02.root 
hadd files_initial/ZJetsToNuNu_HT2500-Inf_final.root rootFiles/ZJetsToNuNu_HT2500-Inf_00.root 
hadd files_initial/ZJetsToNuNu_HT400-600_final.root rootFiles/ZJetsToNuNu_HT400-600_00.root 
hadd files_initial/ZJetsToNuNu_HT800-1200_final.root rootFiles/ZJetsToNuNu_HT800-1200_00.root 
hadd files_initial/ZZTo2L2Q_final.root rootFiles/ZZTo2L2Q_00.root rootFiles/ZZTo2L2Q_01.root rootFiles/ZZTo2L2Q_02.root 


