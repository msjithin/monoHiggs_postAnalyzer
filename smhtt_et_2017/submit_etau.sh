
outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

###########################   MC  #########################

./rootcom etau_analyzer executable_etau  


./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_ext1_v1_00.root DYJetsToLL_M-50_TuneCP5_ext1_v1_00.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_ext1_v1_01.root DYJetsToLL_M-50_TuneCP5_ext1_v1_01.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_ext1_v1_02.root DYJetsToLL_M-50_TuneCP5_ext1_v1_02.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_ext1_v1_03.root DYJetsToLL_M-50_TuneCP5_ext1_v1_03.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_ext1_v1_04.root DYJetsToLL_M-50_TuneCP5_ext1_v1_04.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_ext1_v1_05.root DYJetsToLL_M-50_TuneCP5_ext1_v1_05.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_ext1_v1_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_v1_00.root DYJetsToLL_M-50_TuneCP5_v1_00.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_v1_01.root DYJetsToLL_M-50_TuneCP5_v1_01.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_v1_02.root DYJetsToLL_M-50_TuneCP5_v1_02.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_v1_03.root DYJetsToLL_M-50_TuneCP5_v1_03.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_v1_04.root DYJetsToLL_M-50_TuneCP5_v1_04.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/DYJetsToLL_M-50_TuneCP5_v1_05.root DYJetsToLL_M-50_TuneCP5_v1_05.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY1JetsToLL_M-50_TuneCP5_00.root DY1JetsToLL_M-50_TuneCP5_00.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY1JetsToLL_M-50_TuneCP5_01.root DY1JetsToLL_M-50_TuneCP5_01.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY1JetsToLL_M-50_TuneCP5_02.root DY1JetsToLL_M-50_TuneCP5_02.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY1JetsToLL_M-50_TuneCP5_03.root DY1JetsToLL_M-50_TuneCP5_03.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY1JetsToLL_M-50_TuneCP5_04.root DY1JetsToLL_M-50_TuneCP5_04.root -1 1000 2017 MC DY1JetsToLL_M-50_TuneCP5_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY2JetsToLL_M-50_TuneCP5_00.root DY2JetsToLL_M-50_TuneCP5_00.root -1 1000 2017 MC DY2JetsToLL_M-50_TuneCP5_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY2JetsToLL_M-50_TuneCP5_01.root DY2JetsToLL_M-50_TuneCP5_01.root -1 1000 2017 MC DY2JetsToLL_M-50_TuneCP5_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY3JetsToLL_M-50_TuneCP5_ext1_00.root DY3JetsToLL_M-50_TuneCP5_ext1_00.root -1 1000 2017 MC DY3JetsToLL_M-50_TuneCP5_ext1_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY3JetsToLL_M-50_TuneCP5_v1_00.root DY3JetsToLL_M-50_TuneCP5_v1_00.root -1 1000 2017 MC DY3JetsToLL_M-50_TuneCP5_v1_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/DY4JetsToLL_M-50_TuneCP5_00.root DY4JetsToLL_M-50_TuneCP5_00.root -1 1000 2017 MC DY4JetsToLL_M-50_TuneCP5_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTTo2L2Nu_TuneCP5_00.root TTTo2L2Nu_TuneCP5_00.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTTo2L2Nu_TuneCP5_01.root TTTo2L2Nu_TuneCP5_01.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTTo2L2Nu_TuneCP5_02.root TTTo2L2Nu_TuneCP5_02.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTTo2L2Nu_TuneCP5_03.root TTTo2L2Nu_TuneCP5_03.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTTo2L2Nu_TuneCP5_04.root TTTo2L2Nu_TuneCP5_04.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTTo2L2Nu_TuneCP5_05.root TTTo2L2Nu_TuneCP5_05.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_05 $outDir


./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToHadronic_TuneCP5_00.root TTToHadronic_TuneCP5_00.root -1 1000 2017 MC TTToHadronic_TuneCP5_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToHadronic_TuneCP5_01.root TTToHadronic_TuneCP5_01.root -1 1000 2017 MC TTToHadronic_TuneCP5_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToHadronic_TuneCP5_02.root TTToHadronic_TuneCP5_02.root -1 1000 2017 MC TTToHadronic_TuneCP5_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToHadronic_TuneCP5_03.root TTToHadronic_TuneCP5_03.root -1 1000 2017 MC TTToHadronic_TuneCP5_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToHadronic_TuneCP5_04.root TTToHadronic_TuneCP5_04.root -1 1000 2017 MC TTToHadronic_TuneCP5_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToHadronic_TuneCP5_05.root TTToHadronic_TuneCP5_05.root -1 1000 2017 MC TTToHadronic_TuneCP5_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToHadronic_TuneCP5_06.root TTToHadronic_TuneCP5_06.root -1 1000 2017 MC TTToHadronic_TuneCP5_06 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToHadronic_TuneCP5_07.root TTToHadronic_TuneCP5_07.root -1 1000 2017 MC TTToHadronic_TuneCP5_07 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToHadronic_TuneCP5_08.root TTToHadronic_TuneCP5_08.root -1 1000 2017 MC TTToHadronic_TuneCP5_08 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToSemiLeptonic_TuneCP5_00.root TTToSemiLeptonic_TuneCP5_00.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToSemiLeptonic_TuneCP5_01.root TTToSemiLeptonic_TuneCP5_01.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToSemiLeptonic_TuneCP5_02.root TTToSemiLeptonic_TuneCP5_02.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToSemiLeptonic_TuneCP5_03.root TTToSemiLeptonic_TuneCP5_03.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/TTToSemiLeptonic_TuneCP5_04.root TTToSemiLeptonic_TuneCP5_04.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_04 $outDir



./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraB_00.root SingleElectron_EraB_00.root -1 1000 2017 DATA SingleElectron_EraB_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraB_01.root SingleElectron_EraB_01.root -1 1000 2017 DATA SingleElectron_EraB_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraC_00.root SingleElectron_EraC_00.root -1 1000 2017 DATA SingleElectron_EraC_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraC_01.root SingleElectron_EraC_01.root -1 1000 2017 DATA SingleElectron_EraC_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraC_02.root SingleElectron_EraC_02.root -1 1000 2017 DATA SingleElectron_EraC_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraC_03.root SingleElectron_EraC_03.root -1 1000 2017 DATA SingleElectron_EraC_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraC_04.root SingleElectron_EraC_04.root -1 1000 2017 DATA SingleElectron_EraC_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraD_00.root SingleElectron_EraD_00.root -1 1000 2017 DATA SingleElectron_EraD_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraD_01.root SingleElectron_EraD_01.root -1 1000 2017 DATA SingleElectron_EraD_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraD_02.root SingleElectron_EraD_02.root -1 1000 2017 DATA SingleElectron_EraD_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraE_00.root SingleElectron_EraE_00.root -1 1000 2017 DATA SingleElectron_EraE_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraE_01.root SingleElectron_EraE_01.root -1 1000 2017 DATA SingleElectron_EraE_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraE_02.root SingleElectron_EraE_02.root -1 1000 2017 DATA SingleElectron_EraE_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraE_03.root SingleElectron_EraE_03.root -1 1000 2017 DATA SingleElectron_EraE_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraE_04.root SingleElectron_EraE_04.root -1 1000 2017 DATA SingleElectron_EraE_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraF_00.root SingleElectron_EraF_00.root -1 1000 2017 DATA SingleElectron_EraF_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraF_01.root SingleElectron_EraF_01.root -1 1000 2017 DATA SingleElectron_EraF_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraF_02.root SingleElectron_EraF_02.root -1 1000 2017 DATA SingleElectron_EraF_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraF_03.root SingleElectron_EraF_03.root -1 1000 2017 DATA SingleElectron_EraF_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraF_04.root SingleElectron_EraF_04.root -1 1000 2017 DATA SingleElectron_EraF_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2017_skimmed/etau/SingleElectron_EraF_05.root SingleElectron_EraF_05.root -1 1000 2017 DATA SingleElectron_EraF_05 $outDir
