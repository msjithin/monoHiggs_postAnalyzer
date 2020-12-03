
outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

###########################   MC  #########################

./rootcom tautau_analyzer analyze_tautau  


./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY1JetsToLL_M-50_TuneCP5_00.root DY1JetsToLL_M-50_TuneCP5_00.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY1JetsToLL_M-50_TuneCP5_01.root DY1JetsToLL_M-50_TuneCP5_01.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY1JetsToLL_M-50_TuneCP5_02.root DY1JetsToLL_M-50_TuneCP5_02.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY1JetsToLL_M-50_TuneCP5_03.root DY1JetsToLL_M-50_TuneCP5_03.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY1JetsToLL_M-50_TuneCP5_04.root DY1JetsToLL_M-50_TuneCP5_04.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_04 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY2JetsToLL_M-50_TuneCP5_00.root DY2JetsToLL_M-50_TuneCP5_00.root -1 1000 2017 MC DY2JetsToLL_M-50_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY2JetsToLL_M-50_TuneCP5_01.root DY2JetsToLL_M-50_TuneCP5_01.root -1 1000 2017 MC DY2JetsToLL_M-50_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY3JetsToLL_M-50_TuneCP5_ext1_00.root DY3JetsToLL_M-50_TuneCP5_ext1_00.root -1 1000 2017 MC DY3JetsToLL_M-50_TuneCP5_ext1_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY3JetsToLL_M-50_TuneCP5_v1_00.root DY3JetsToLL_M-50_TuneCP5_v1_00.root -1 1000 2017 MC DY3JetsToLL_M-50_TuneCP5_v1_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DY4JetsToLL_M-50_TuneCP5_00.root DY4JetsToLL_M-50_TuneCP5_00.root -1 1000 2017 MC DY4JetsToLL_M-50_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-10to50_TuneCP5_00.root DYJetsToLL_M-10to50_TuneCP5_00.root -1 1000 2017 MC DYJetsToLL_M-10to50_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-10to50_TuneCP5_01.root DYJetsToLL_M-10to50_TuneCP5_01.root -1 1000 2017 MC DYJetsToLL_M-10to50_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-10to50_TuneCP5_02.root DYJetsToLL_M-10to50_TuneCP5_02.root -1 1000 2017 MC DYJetsToLL_M-10to50_TuneCP5_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-10to50_TuneCP5_03.root DYJetsToLL_M-10to50_TuneCP5_03.root -1 1000 2017 MC DYJetsToLL_M-10to50_TuneCP5_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-10to50_TuneCP5_04.root DYJetsToLL_M-10to50_TuneCP5_04.root -1 1000 2017 MC DYJetsToLL_M-10to50_TuneCP5_04 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_ext1_v1_00.root DYJetsToLL_M-50_TuneCP5_ext1_v1_00.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_ext1_v1_01.root DYJetsToLL_M-50_TuneCP5_ext1_v1_01.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_ext1_v1_02.root DYJetsToLL_M-50_TuneCP5_ext1_v1_02.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_ext1_v1_03.root DYJetsToLL_M-50_TuneCP5_ext1_v1_03.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_ext1_v1_04.root DYJetsToLL_M-50_TuneCP5_ext1_v1_04.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_04 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_ext1_v1_05.root DYJetsToLL_M-50_TuneCP5_ext1_v1_05.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_05 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_v1_00.root DYJetsToLL_M-50_TuneCP5_v1_00.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_v1_01.root DYJetsToLL_M-50_TuneCP5_v1_01.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_v1_02.root DYJetsToLL_M-50_TuneCP5_v1_02.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_v1_03.root DYJetsToLL_M-50_TuneCP5_v1_03.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_v1_04.root DYJetsToLL_M-50_TuneCP5_v1_04.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_04 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/DYJetsToLL_M-50_TuneCP5_v1_05.root DYJetsToLL_M-50_TuneCP5_v1_05.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_05 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_00.root EWKWMinus2Jets_WToLNu_M-50_TuneCP5_00.root -1 1000 2017 MC EWKWMinus2Jets_WToLNu_M-50_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/EWKWPlus2Jets_WToLNu_M-50_TuneCP5_00.root EWKWPlus2Jets_WToLNu_M-50_TuneCP5_00.root -1 1000 2017 MC EWKWPlus2Jets_WToLNu_M-50_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/EWKZ2Jets_ZToLL_M-50_TuneCP5_00.root EWKZ2Jets_ZToLL_M-50_TuneCP5_00.root -1 1000 2017 MC EWKZ2Jets_ZToLL_M-50_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/EWKZ2Jets_ZToNuNu_TuneCP5_00.root EWKZ2Jets_ZToNuNu_TuneCP5_00.root -1 1000 2017 MC EWKZ2Jets_ZToNuNu_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/GluGluHToTauTau_M125_00.root GluGluHToTauTau_M125_00.root -1 1000 2017 MC GluGluHToTauTau_M125_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/GluGluHToWWTo2L2Nu_M125_00.root GluGluHToWWTo2L2Nu_M125_00.root -1 1000 2017 MC GluGluHToWWTo2L2Nu_M125_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_00.root ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_00.root -1 1000 2017 MC ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_00.root ST_t-channel_top_4f_inclusiveDecays_TuneCP5_00.root -1 1000 2017 MC ST_t-channel_top_4f_inclusiveDecays_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_00.root ST_tW_antitop_5f_inclusiveDecays_TuneCP5_00.root -1 1000 2017 MC ST_tW_antitop_5f_inclusiveDecays_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ST_tW_top_5f_inclusiveDecays_TuneCP5_00.root ST_tW_top_5f_inclusiveDecays_TuneCP5_00.root -1 1000 2017 MC ST_tW_top_5f_inclusiveDecays_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTTo2L2Nu_TuneCP5_00.root TTTo2L2Nu_TuneCP5_00.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTTo2L2Nu_TuneCP5_01.root TTTo2L2Nu_TuneCP5_01.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTTo2L2Nu_TuneCP5_02.root TTTo2L2Nu_TuneCP5_02.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTTo2L2Nu_TuneCP5_03.root TTTo2L2Nu_TuneCP5_03.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTTo2L2Nu_TuneCP5_04.root TTTo2L2Nu_TuneCP5_04.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_04 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTTo2L2Nu_TuneCP5_05.root TTTo2L2Nu_TuneCP5_05.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_05 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToHadronic_TuneCP5_00.root TTToHadronic_TuneCP5_00.root -1 1000 2017 MC TTToHadronic_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToHadronic_TuneCP5_01.root TTToHadronic_TuneCP5_01.root -1 1000 2017 MC TTToHadronic_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToHadronic_TuneCP5_02.root TTToHadronic_TuneCP5_02.root -1 1000 2017 MC TTToHadronic_TuneCP5_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToHadronic_TuneCP5_03.root TTToHadronic_TuneCP5_03.root -1 1000 2017 MC TTToHadronic_TuneCP5_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToHadronic_TuneCP5_04.root TTToHadronic_TuneCP5_04.root -1 1000 2017 MC TTToHadronic_TuneCP5_04 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToHadronic_TuneCP5_05.root TTToHadronic_TuneCP5_05.root -1 1000 2017 MC TTToHadronic_TuneCP5_05 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToHadronic_TuneCP5_06.root TTToHadronic_TuneCP5_06.root -1 1000 2017 MC TTToHadronic_TuneCP5_06 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToHadronic_TuneCP5_07.root TTToHadronic_TuneCP5_07.root -1 1000 2017 MC TTToHadronic_TuneCP5_07 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToHadronic_TuneCP5_08.root TTToHadronic_TuneCP5_08.root -1 1000 2017 MC TTToHadronic_TuneCP5_08 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToSemiLeptonic_TuneCP5_00.root TTToSemiLeptonic_TuneCP5_00.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToSemiLeptonic_TuneCP5_01.root TTToSemiLeptonic_TuneCP5_01.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToSemiLeptonic_TuneCP5_02.root TTToSemiLeptonic_TuneCP5_02.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToSemiLeptonic_TuneCP5_03.root TTToSemiLeptonic_TuneCP5_03.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToSemiLeptonic_TuneCP5_04.root TTToSemiLeptonic_TuneCP5_04.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_04 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToSemiLeptonic_TuneCP5_05.root TTToSemiLeptonic_TuneCP5_05.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_05 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToSemiLeptonic_TuneCP5_06.root TTToSemiLeptonic_TuneCP5_06.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_06 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToSemiLeptonic_TuneCP5_07.root TTToSemiLeptonic_TuneCP5_07.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_07 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/TTToSemiLeptonic_TuneCP5_08.root TTToSemiLeptonic_TuneCP5_08.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_08 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/VBFHToTauTau_M125_00.root VBFHToTauTau_M125_00.root -1 1000 2017 MC VBFHToTauTau_M125_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/VBFHToWWTo2L2Nu_M125_00.root VBFHToWWTo2L2Nu_M125_00.root -1 1000 2017 MC VBFHToWWTo2L2Nu_M125_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/VVTo2L2Nu_00.root VVTo2L2Nu_00.root -1 1000 2017 MC VVTo2L2Nu_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W1JetsToLNu_TuneCP5_00.root W1JetsToLNu_TuneCP5_00.root -1 1000 2017 MC W1JetsToLNu_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W1JetsToLNu_TuneCP5_01.root W1JetsToLNu_TuneCP5_01.root -1 1000 2017 MC W1JetsToLNu_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W1JetsToLNu_TuneCP5_02.root W1JetsToLNu_TuneCP5_02.root -1 1000 2017 MC W1JetsToLNu_TuneCP5_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W1JetsToLNu_TuneCP5_03.root W1JetsToLNu_TuneCP5_03.root -1 1000 2017 MC W1JetsToLNu_TuneCP5_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W1JetsToLNu_TuneCP5_04.root W1JetsToLNu_TuneCP5_04.root -1 1000 2017 MC W1JetsToLNu_TuneCP5_04 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W1JetsToLNu_TuneCP5_05.root W1JetsToLNu_TuneCP5_05.root -1 1000 2017 MC W1JetsToLNu_TuneCP5_05 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W1JetsToLNu_TuneCP5_06.root W1JetsToLNu_TuneCP5_06.root -1 1000 2017 MC W1JetsToLNu_TuneCP5_06 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W2JetsToLNu_TuneCP5_00.root W2JetsToLNu_TuneCP5_00.root -1 1000 2017 MC W2JetsToLNu_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W3JetsToLNu_TuneCP5_00.root W3JetsToLNu_TuneCP5_00.root -1 1000 2017 MC W3JetsToLNu_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W3JetsToLNu_TuneCP5_01.root W3JetsToLNu_TuneCP5_01.root -1 1000 2017 MC W3JetsToLNu_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W4JetsToLNu_TuneCP5_00.root W4JetsToLNu_TuneCP5_00.root -1 1000 2017 MC W4JetsToLNu_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/W4JetsToLNu_TuneCP5_01.root W4JetsToLNu_TuneCP5_01.root -1 1000 2017 MC W4JetsToLNu_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WJetsToLNu_TuneCP5_00.root WJetsToLNu_TuneCP5_00.root -1 1000 2017 MC WJetsToLNu_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WJetsToLNu_TuneCP5_01.root WJetsToLNu_TuneCP5_01.root -1 1000 2017 MC WJetsToLNu_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WJetsToLNu_TuneCP5_02.root WJetsToLNu_TuneCP5_02.root -1 1000 2017 MC WJetsToLNu_TuneCP5_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WJetsToLNu_TuneCP5_03.root WJetsToLNu_TuneCP5_03.root -1 1000 2017 MC WJetsToLNu_TuneCP5_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WJetsToLNu_TuneCP5_04.root WJetsToLNu_TuneCP5_04.root -1 1000 2017 MC WJetsToLNu_TuneCP5_04 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WWTo1L1Nu2Q_00.root WWTo1L1Nu2Q_00.root -1 1000 2017 MC WWTo1L1Nu2Q_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WWToLNuQQ_NNPDF31_TuneCP5_00.root WWToLNuQQ_NNPDF31_TuneCP5_00.root -1 1000 2017 MC WWToLNuQQ_NNPDF31_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WWW_4F_TuneCP5_00.root WWW_4F_TuneCP5_00.root -1 1000 2017 MC WWW_4F_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WWZ_4F_TuneCP5_00.root WWZ_4F_TuneCP5_00.root -1 1000 2017 MC WWZ_4F_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WW_TuneCP5_00.root WW_TuneCP5_00.root -1 1000 2017 MC WW_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WZTo3LNu_TuneCP5_00.root WZTo3LNu_TuneCP5_00.root -1 1000 2017 MC WZTo3LNu_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WZTo3LNu_TuneCP5_01.root WZTo3LNu_TuneCP5_01.root -1 1000 2017 MC WZTo3LNu_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WZZ_TuneCP5_00.root WZZ_TuneCP5_00.root -1 1000 2017 MC WZZ_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WZ_TuneCP5_00.root WZ_TuneCP5_00.root -1 1000 2017 MC WZ_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WminusHToTauTau_M125_00.root WminusHToTauTau_M125_00.root -1 1000 2017 MC WminusHToTauTau_M125_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/WplusHToTauTau_M125_00.root WplusHToTauTau_M125_00.root -1 1000 2017 MC WplusHToTauTau_M125_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZHToTauTau_M125_00.root ZHToTauTau_M125_00.root -1 1000 2017 MC ZHToTauTau_M125_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-100To200_00.root ZJetsToNuNu_HT-100To200_00.root -1 1000 2017 MC ZJetsToNuNu_HT-100To200_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-100To200_01.root ZJetsToNuNu_HT-100To200_01.root -1 1000 2017 MC ZJetsToNuNu_HT-100To200_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-100To200_02.root ZJetsToNuNu_HT-100To200_02.root -1 1000 2017 MC ZJetsToNuNu_HT-100To200_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-1200To2500_00.root ZJetsToNuNu_HT-1200To2500_00.root -1 1000 2017 MC ZJetsToNuNu_HT-1200To2500_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-200To400_00.root ZJetsToNuNu_HT-200To400_00.root -1 1000 2017 MC ZJetsToNuNu_HT-200To400_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-200To400_01.root ZJetsToNuNu_HT-200To400_01.root -1 1000 2017 MC ZJetsToNuNu_HT-200To400_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-200To400_02.root ZJetsToNuNu_HT-200To400_02.root -1 1000 2017 MC ZJetsToNuNu_HT-200To400_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-2500ToInf_00.root ZJetsToNuNu_HT-2500ToInf_00.root -1 1000 2017 MC ZJetsToNuNu_HT-2500ToInf_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-400To600_00.root ZJetsToNuNu_HT-400To600_00.root -1 1000 2017 MC ZJetsToNuNu_HT-400To600_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-600To800_00.root ZJetsToNuNu_HT-600To800_00.root -1 1000 2017 MC ZJetsToNuNu_HT-600To800_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZJetsToNuNu_HT-800To1200_00.root ZJetsToNuNu_HT-800To1200_00.root -1 1000 2017 MC ZJetsToNuNu_HT-800To1200_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZZTo2L2Q_00.root ZZTo2L2Q_00.root -1 1000 2017 MC ZZTo2L2Q_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZZTo2L2Q_01.root ZZTo2L2Q_01.root -1 1000 2017 MC ZZTo2L2Q_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZZTo2L2Q_02.root ZZTo2L2Q_02.root -1 1000 2017 MC ZZTo2L2Q_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZZTo4L_TuneCP5_00.root ZZTo4L_TuneCP5_00.root -1 1000 2017 MC ZZTo4L_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZZTo4L_TuneCP5_01.root ZZTo4L_TuneCP5_01.root -1 1000 2017 MC ZZTo4L_TuneCP5_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZZTo4L_TuneCP5_02.root ZZTo4L_TuneCP5_02.root -1 1000 2017 MC ZZTo4L_TuneCP5_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZZZ_TuneCP5_00.root ZZZ_TuneCP5_00.root -1 1000 2017 MC ZZZ_TuneCP5_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/ZZ_TuneCP5_00.root ZZ_TuneCP5_00.root -1 1000 2017 MC ZZ_TuneCP5_00 $outDir


###########################  DATA #########################

./rootcom tautau_analyzer analyze_tautau 


./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraB_00.root Tau_EraB_00.root -1 1000 2017 DATA Tau_EraB_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraB_01.root Tau_EraB_01.root -1 1000 2017 DATA Tau_EraB_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraC_00.root Tau_EraC_00.root -1 1000 2017 DATA Tau_EraC_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraC_01.root Tau_EraC_01.root -1 1000 2017 DATA Tau_EraC_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraC_02.root Tau_EraC_02.root -1 1000 2017 DATA Tau_EraC_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraC_03.root Tau_EraC_03.root -1 1000 2017 DATA Tau_EraC_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraD_00.root Tau_EraD_00.root -1 1000 2017 DATA Tau_EraD_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraD_01.root Tau_EraD_01.root -1 1000 2017 DATA Tau_EraD_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraD_02.root Tau_EraD_02.root -1 1000 2017 DATA Tau_EraD_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraE_00.root Tau_EraE_00.root -1 1000 2017 DATA Tau_EraE_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraE_01.root Tau_EraE_01.root -1 1000 2017 DATA Tau_EraE_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraE_02.root Tau_EraE_02.root -1 1000 2017 DATA Tau_EraE_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraE_03.root Tau_EraE_03.root -1 1000 2017 DATA Tau_EraE_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraF_00.root Tau_EraF_00.root -1 1000 2017 DATA Tau_EraF_00 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraF_01.root Tau_EraF_01.root -1 1000 2017 DATA Tau_EraF_01 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraF_02.root Tau_EraF_02.root -1 1000 2017 DATA Tau_EraF_02 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraF_03.root Tau_EraF_03.root -1 1000 2017 DATA Tau_EraF_03 $outDir
./MakeCondorFiles.csh analyze_tautau /hdfs/store/user/jmadhusu/2017_skimmed/tautau_new/Tau_EraF_04.root Tau_EraF_04.root -1 1000 2017 DATA Tau_EraF_04 $outDir
