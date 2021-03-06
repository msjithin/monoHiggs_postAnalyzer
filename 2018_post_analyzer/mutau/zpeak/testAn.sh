
outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

###########################   MC  #########################

./rootcom mutau_analyzer analyze_mutau  


./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY1JetsToLL_00.root DY1JetsToLL_00.root -1 1000 2018_test MC DY1JetsToLL_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY1JetsToLL_01.root DY1JetsToLL_01.root -1 1000 2018_test MC DY1JetsToLL_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY1JetsToLL_02.root DY1JetsToLL_02.root -1 1000 2018_test MC DY1JetsToLL_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY1JetsToLL_03.root DY1JetsToLL_03.root -1 1000 2018_test MC DY1JetsToLL_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY1JetsToLL_04.root DY1JetsToLL_04.root -1 1000 2018_test MC DY1JetsToLL_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY1JetsToLL_05.root DY1JetsToLL_05.root -1 1000 2018_test MC DY1JetsToLL_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY2JetsToLL_00.root DY2JetsToLL_00.root -1 1000 2018_test MC DY2JetsToLL_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY2JetsToLL_01.root DY2JetsToLL_01.root -1 1000 2018_test MC DY2JetsToLL_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY2JetsToLL_02.root DY2JetsToLL_02.root -1 1000 2018_test MC DY2JetsToLL_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY3JetsToLL_00.root DY3JetsToLL_00.root -1 1000 2018_test MC DY3JetsToLL_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DY4JetsToLL_00.root DY4JetsToLL_00.root -1 1000 2018_test MC DY4JetsToLL_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_00.root DYJetsToLL_00.root -1 1000 2018_test MC DYJetsToLL_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_01.root DYJetsToLL_01.root -1 1000 2018_test MC DYJetsToLL_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_02.root DYJetsToLL_02.root -1 1000 2018_test MC DYJetsToLL_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_03.root DYJetsToLL_03.root -1 1000 2018_test MC DYJetsToLL_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_04.root DYJetsToLL_04.root -1 1000 2018_test MC DYJetsToLL_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_05.root DYJetsToLL_05.root -1 1000 2018_test MC DYJetsToLL_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_06.root DYJetsToLL_06.root -1 1000 2018_test MC DYJetsToLL_06 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_M10to50_00.root DYJetsToLL_M10to50_00.root -1 1000 2018_test MC DYJetsToLL_M10to50_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_M10to50_01.root DYJetsToLL_M10to50_01.root -1 1000 2018_test MC DYJetsToLL_M10to50_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_M10to50_02.root DYJetsToLL_M10to50_02.root -1 1000 2018_test MC DYJetsToLL_M10to50_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/DYJetsToLL_M10to50_03.root DYJetsToLL_M10to50_03.root -1 1000 2018_test MC DYJetsToLL_M10to50_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/EWKWMinus2Jets_WToLNu_00.root EWKWMinus2Jets_WToLNu_00.root -1 1000 2018_test MC EWKWMinus2Jets_WToLNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/EWKWPlus2Jets_WToLNu_00.root EWKWPlus2Jets_WToLNu_00.root -1 1000 2018_test MC EWKWPlus2Jets_WToLNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/EWKZ2Jets_ZToLL_00.root EWKZ2Jets_ZToLL_00.root -1 1000 2018_test MC EWKZ2Jets_ZToLL_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/EWKZ2Jets_ZToNuNu_00.root EWKZ2Jets_ZToNuNu_00.root -1 1000 2018_test MC EWKZ2Jets_ZToNuNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/GluGluHToTauTau_00.root GluGluHToTauTau_00.root -1 1000 2018_test MC GluGluHToTauTau_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/GluGluHToTauTau_01.root GluGluHToTauTau_01.root -1 1000 2018_test MC GluGluHToTauTau_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/GluGluHToWWTo2L2Nu_00.root GluGluHToWWTo2L2Nu_00.root -1 1000 2018_test MC GluGluHToWWTo2L2Nu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/GluGluZH_HToWW_00.root GluGluZH_HToWW_00.root -1 1000 2018_test MC GluGluZH_HToWW_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/HWminusJ_HToWW_00.root HWminusJ_HToWW_00.root -1 1000 2018_test MC HWminusJ_HToWW_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/HWplusJ_HToWW_00.root HWplusJ_HToWW_00.root -1 1000 2018_test MC HWplusJ_HToWW_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/HZJ_HToWW_00.root HZJ_HToWW_00.root -1 1000 2018_test MC HZJ_HToWW_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_antitop_00.root ST_t-channel_antitop_00.root -1 1000 2018_test MC ST_t-channel_antitop_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_antitop_01.root ST_t-channel_antitop_01.root -1 1000 2018_test MC ST_t-channel_antitop_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_antitop_02.root ST_t-channel_antitop_02.root -1 1000 2018_test MC ST_t-channel_antitop_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_antitop_03.root ST_t-channel_antitop_03.root -1 1000 2018_test MC ST_t-channel_antitop_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_antitop_04.root ST_t-channel_antitop_04.root -1 1000 2018_test MC ST_t-channel_antitop_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_antitop_05.root ST_t-channel_antitop_05.root -1 1000 2018_test MC ST_t-channel_antitop_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_antitop_06.root ST_t-channel_antitop_06.root -1 1000 2018_test MC ST_t-channel_antitop_06 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_antitop_07.root ST_t-channel_antitop_07.root -1 1000 2018_test MC ST_t-channel_antitop_07 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_00.root ST_t-channel_top_00.root -1 1000 2018_test MC ST_t-channel_top_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_01.root ST_t-channel_top_01.root -1 1000 2018_test MC ST_t-channel_top_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_02.root ST_t-channel_top_02.root -1 1000 2018_test MC ST_t-channel_top_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_03.root ST_t-channel_top_03.root -1 1000 2018_test MC ST_t-channel_top_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_04.root ST_t-channel_top_04.root -1 1000 2018_test MC ST_t-channel_top_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_05.root ST_t-channel_top_05.root -1 1000 2018_test MC ST_t-channel_top_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_06.root ST_t-channel_top_06.root -1 1000 2018_test MC ST_t-channel_top_06 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_07.root ST_t-channel_top_07.root -1 1000 2018_test MC ST_t-channel_top_07 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_08.root ST_t-channel_top_08.root -1 1000 2018_test MC ST_t-channel_top_08 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_09.root ST_t-channel_top_09.root -1 1000 2018_test MC ST_t-channel_top_09 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_t-channel_top_10.root ST_t-channel_top_10.root -1 1000 2018_test MC ST_t-channel_top_10 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_tW_antitop_00.root ST_tW_antitop_00.root -1 1000 2018_test MC ST_tW_antitop_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ST_tW_top_00.root ST_tW_top_00.root -1 1000 2018_test MC ST_tW_top_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTTo2L2Nu_powheg_00.root TTTo2L2Nu_powheg_00.root -1 1000 2018_test MC TTTo2L2Nu_powheg_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTTo2L2Nu_powheg_01.root TTTo2L2Nu_powheg_01.root -1 1000 2018_test MC TTTo2L2Nu_powheg_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTTo2L2Nu_powheg_02.root TTTo2L2Nu_powheg_02.root -1 1000 2018_test MC TTTo2L2Nu_powheg_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTTo2L2Nu_powheg_03.root TTTo2L2Nu_powheg_03.root -1 1000 2018_test MC TTTo2L2Nu_powheg_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTTo2L2Nu_powheg_04.root TTTo2L2Nu_powheg_04.root -1 1000 2018_test MC TTTo2L2Nu_powheg_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTTo2L2Nu_powheg_05.root TTTo2L2Nu_powheg_05.root -1 1000 2018_test MC TTTo2L2Nu_powheg_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTTo2L2Nu_powheg_06.root TTTo2L2Nu_powheg_06.root -1 1000 2018_test MC TTTo2L2Nu_powheg_06 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_00.root TTToHadronic_powheg_00.root -1 1000 2018_test MC TTToHadronic_powheg_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_01.root TTToHadronic_powheg_01.root -1 1000 2018_test MC TTToHadronic_powheg_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_02.root TTToHadronic_powheg_02.root -1 1000 2018_test MC TTToHadronic_powheg_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_03.root TTToHadronic_powheg_03.root -1 1000 2018_test MC TTToHadronic_powheg_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_04.root TTToHadronic_powheg_04.root -1 1000 2018_test MC TTToHadronic_powheg_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_05.root TTToHadronic_powheg_05.root -1 1000 2018_test MC TTToHadronic_powheg_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_06.root TTToHadronic_powheg_06.root -1 1000 2018_test MC TTToHadronic_powheg_06 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_07.root TTToHadronic_powheg_07.root -1 1000 2018_test MC TTToHadronic_powheg_07 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_08.root TTToHadronic_powheg_08.root -1 1000 2018_test MC TTToHadronic_powheg_08 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_09.root TTToHadronic_powheg_09.root -1 1000 2018_test MC TTToHadronic_powheg_09 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_10.root TTToHadronic_powheg_10.root -1 1000 2018_test MC TTToHadronic_powheg_10 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToHadronic_powheg_11.root TTToHadronic_powheg_11.root -1 1000 2018_test MC TTToHadronic_powheg_11 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToSemiLeptonic_powheg_00.root TTToSemiLeptonic_powheg_00.root -1 1000 2018_test MC TTToSemiLeptonic_powheg_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToSemiLeptonic_powheg_01.root TTToSemiLeptonic_powheg_01.root -1 1000 2018_test MC TTToSemiLeptonic_powheg_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToSemiLeptonic_powheg_02.root TTToSemiLeptonic_powheg_02.root -1 1000 2018_test MC TTToSemiLeptonic_powheg_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToSemiLeptonic_powheg_03.root TTToSemiLeptonic_powheg_03.root -1 1000 2018_test MC TTToSemiLeptonic_powheg_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToSemiLeptonic_powheg_04.root TTToSemiLeptonic_powheg_04.root -1 1000 2018_test MC TTToSemiLeptonic_powheg_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToSemiLeptonic_powheg_05.root TTToSemiLeptonic_powheg_05.root -1 1000 2018_test MC TTToSemiLeptonic_powheg_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/TTToSemiLeptonic_powheg_06.root TTToSemiLeptonic_powheg_06.root -1 1000 2018_test MC TTToSemiLeptonic_powheg_06 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/VBFHToTauTau_00.root VBFHToTauTau_00.root -1 1000 2018_test MC VBFHToTauTau_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/VBFHToWWTo2L2Nu_00.root VBFHToWWTo2L2Nu_00.root -1 1000 2018_test MC VBFHToWWTo2L2Nu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/VVTo2L2Nu_00.root VVTo2L2Nu_00.root -1 1000 2018_test MC VVTo2L2Nu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/VVTo2L2Nu_01.root VVTo2L2Nu_01.root -1 1000 2018_test MC VVTo2L2Nu_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W1JetsToLNu_00.root W1JetsToLNu_00.root -1 1000 2018_test MC W1JetsToLNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W1JetsToLNu_01.root W1JetsToLNu_01.root -1 1000 2018_test MC W1JetsToLNu_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W1JetsToLNu_02.root W1JetsToLNu_02.root -1 1000 2018_test MC W1JetsToLNu_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W1JetsToLNu_03.root W1JetsToLNu_03.root -1 1000 2018_test MC W1JetsToLNu_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W1JetsToLNu_04.root W1JetsToLNu_04.root -1 1000 2018_test MC W1JetsToLNu_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W1JetsToLNu_05.root W1JetsToLNu_05.root -1 1000 2018_test MC W1JetsToLNu_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W2JetsToLNu_00.root W2JetsToLNu_00.root -1 1000 2018_test MC W2JetsToLNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W2JetsToLNu_01.root W2JetsToLNu_01.root -1 1000 2018_test MC W2JetsToLNu_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W2JetsToLNu_02.root W2JetsToLNu_02.root -1 1000 2018_test MC W2JetsToLNu_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W3JetsToLNu_00.root W3JetsToLNu_00.root -1 1000 2018_test MC W3JetsToLNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W3JetsToLNu_01.root W3JetsToLNu_01.root -1 1000 2018_test MC W3JetsToLNu_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W4JetsToLNu_00.root W4JetsToLNu_00.root -1 1000 2018_test MC W4JetsToLNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/W4JetsToLNu_01.root W4JetsToLNu_01.root -1 1000 2018_test MC W4JetsToLNu_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WGToLNuG_00.root WGToLNuG_00.root -1 1000 2018_test MC WGToLNuG_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WJetsToLNu_00.root WJetsToLNu_00.root -1 1000 2018_test MC WJetsToLNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WJetsToLNu_01.root WJetsToLNu_01.root -1 1000 2018_test MC WJetsToLNu_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WJetsToLNu_02.root WJetsToLNu_02.root -1 1000 2018_test MC WJetsToLNu_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WJetsToLNu_03.root WJetsToLNu_03.root -1 1000 2018_test MC WJetsToLNu_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WJetsToLNu_04.root WJetsToLNu_04.root -1 1000 2018_test MC WJetsToLNu_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WJetsToLNu_05.root WJetsToLNu_05.root -1 1000 2018_test MC WJetsToLNu_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WJetsToLNu_06.root WJetsToLNu_06.root -1 1000 2018_test MC WJetsToLNu_06 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WJetsToLNu_07.root WJetsToLNu_07.root -1 1000 2018_test MC WJetsToLNu_07 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WJetsToLNu_08.root WJetsToLNu_08.root -1 1000 2018_test MC WJetsToLNu_08 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WWTo1L1Nu2Q_00.root WWTo1L1Nu2Q_00.root -1 1000 2018_test MC WWTo1L1Nu2Q_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WWToLNuQQ_00.root WWToLNuQQ_00.root -1 1000 2018_test MC WWToLNuQQ_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WWToLNuQQ_01.root WWToLNuQQ_01.root -1 1000 2018_test MC WWToLNuQQ_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WWW_00.root WWW_00.root -1 1000 2018_test MC WWW_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WWZ_00.root WWZ_00.root -1 1000 2018_test MC WWZ_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WZTo3LNu_00.root WZTo3LNu_00.root -1 1000 2018_test MC WZTo3LNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WZTo3LNu_01.root WZTo3LNu_01.root -1 1000 2018_test MC WZTo3LNu_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WZZ_00.root WZZ_00.root -1 1000 2018_test MC WZZ_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WminusHToTauTau_00.root WminusHToTauTau_00.root -1 1000 2018_test MC WminusHToTauTau_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/WplusHToTauTau_00.root WplusHToTauTau_00.root -1 1000 2018_test MC WplusHToTauTau_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZHToTauTau_00.root ZHToTauTau_00.root -1 1000 2018_test MC ZHToTauTau_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT100-200_00.root ZJetsToNuNu_HT100-200_00.root -1 1000 2018_test MC ZJetsToNuNu_HT100-200_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT100-200_01.root ZJetsToNuNu_HT100-200_01.root -1 1000 2018_test MC ZJetsToNuNu_HT100-200_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT100-200_02.root ZJetsToNuNu_HT100-200_02.root -1 1000 2018_test MC ZJetsToNuNu_HT100-200_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT1200-2500_00.root ZJetsToNuNu_HT1200-2500_00.root -1 1000 2018_test MC ZJetsToNuNu_HT1200-2500_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT200-400_00.root ZJetsToNuNu_HT200-400_00.root -1 1000 2018_test MC ZJetsToNuNu_HT200-400_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT200-400_01.root ZJetsToNuNu_HT200-400_01.root -1 1000 2018_test MC ZJetsToNuNu_HT200-400_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT200-400_02.root ZJetsToNuNu_HT200-400_02.root -1 1000 2018_test MC ZJetsToNuNu_HT200-400_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT2500-Inf_00.root ZJetsToNuNu_HT2500-Inf_00.root -1 1000 2018_test MC ZJetsToNuNu_HT2500-Inf_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT400-600_00.root ZJetsToNuNu_HT400-600_00.root -1 1000 2018_test MC ZJetsToNuNu_HT400-600_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT600-800_00.root ZJetsToNuNu_HT600-800_00.root -1 1000 2018_test MC ZJetsToNuNu_HT600-800_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZJetsToNuNu_HT800-1200_00.root ZJetsToNuNu_HT800-1200_00.root -1 1000 2018_test MC ZJetsToNuNu_HT800-1200_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZZTo2L2Q_00.root ZZTo2L2Q_00.root -1 1000 2018_test MC ZZTo2L2Q_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZZTo2L2Q_01.root ZZTo2L2Q_01.root -1 1000 2018_test MC ZZTo2L2Q_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZZTo2L2Q_02.root ZZTo2L2Q_02.root -1 1000 2018_test MC ZZTo2L2Q_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZZTo4L_00.root ZZTo4L_00.root -1 1000 2018_test MC ZZTo4L_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ZZZ_00.root ZZZ_00.root -1 1000 2018_test MC ZZZ_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ggZH_HToTauTau_ZToLL_00.root ggZH_HToTauTau_ZToLL_00.root -1 1000 2018_test MC ggZH_HToTauTau_ZToLL_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ggZH_HToTauTau_ZToNuNu_00.root ggZH_HToTauTau_ZToNuNu_00.root -1 1000 2018_test MC ggZH_HToTauTau_ZToNuNu_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ggZH_HToTauTau_ZToQQ_00.root ggZH_HToTauTau_ZToQQ_00.root -1 1000 2018_test MC ggZH_HToTauTau_ZToQQ_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ttHToNonbb_00.root ttHToNonbb_00.root -1 1000 2018_test MC ttHToNonbb_00 $outDir


###########################  DATA #########################

./rootcom mutau_analyzer analyze_mutau  


./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonA_00.root SingleMuonA_00.root -1 1000 2018_test DATA SingleMuonA_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonA_01.root SingleMuonA_01.root -1 1000 2018_test DATA SingleMuonA_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonA_02.root SingleMuonA_02.root -1 1000 2018_test DATA SingleMuonA_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonA_03.root SingleMuonA_03.root -1 1000 2018_test DATA SingleMuonA_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonA_04.root SingleMuonA_04.root -1 1000 2018_test DATA SingleMuonA_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonB_00.root SingleMuonB_00.root -1 1000 2018_test DATA SingleMuonB_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonB_01.root SingleMuonB_01.root -1 1000 2018_test DATA SingleMuonB_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonB_02.root SingleMuonB_02.root -1 1000 2018_test DATA SingleMuonB_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonC_00.root SingleMuonC_00.root -1 1000 2018_test DATA SingleMuonC_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonC_01.root SingleMuonC_01.root -1 1000 2018_test DATA SingleMuonC_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonC_02.root SingleMuonC_02.root -1 1000 2018_test DATA SingleMuonC_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_00.root SingleMuonD_PromptReco_00.root -1 1000 2018_test DATA SingleMuonD_PromptReco_00 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_01.root SingleMuonD_PromptReco_01.root -1 1000 2018_test DATA SingleMuonD_PromptReco_01 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_02.root SingleMuonD_PromptReco_02.root -1 1000 2018_test DATA SingleMuonD_PromptReco_02 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_03.root SingleMuonD_PromptReco_03.root -1 1000 2018_test DATA SingleMuonD_PromptReco_03 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_04.root SingleMuonD_PromptReco_04.root -1 1000 2018_test DATA SingleMuonD_PromptReco_04 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_05.root SingleMuonD_PromptReco_05.root -1 1000 2018_test DATA SingleMuonD_PromptReco_05 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_06.root SingleMuonD_PromptReco_06.root -1 1000 2018_test DATA SingleMuonD_PromptReco_06 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_07.root SingleMuonD_PromptReco_07.root -1 1000 2018_test DATA SingleMuonD_PromptReco_07 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_08.root SingleMuonD_PromptReco_08.root -1 1000 2018_test DATA SingleMuonD_PromptReco_08 $outDir
./MakeCondorFiles.csh analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/SingleMuonD_PromptReco_09.root SingleMuonD_PromptReco_09.root -1 1000 2018_test DATA SingleMuonD_PromptReco_09 $outDir
