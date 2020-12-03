
outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

###########################   MC  #########################

./rootcom etau_analyzer executable_etau  
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY1JetsToLL_00.root DY1JetsToLL_00.root -1 1000 2018 MC DY1JetsToLL_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY1JetsToLL_01.root DY1JetsToLL_01.root -1 1000 2018 MC DY1JetsToLL_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY1JetsToLL_02.root DY1JetsToLL_02.root -1 1000 2018 MC DY1JetsToLL_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY1JetsToLL_03.root DY1JetsToLL_03.root -1 1000 2018 MC DY1JetsToLL_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY1JetsToLL_04.root DY1JetsToLL_04.root -1 1000 2018 MC DY1JetsToLL_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY2JetsToLL_00.root DY2JetsToLL_00.root -1 1000 2018 MC DY2JetsToLL_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY2JetsToLL_01.root DY2JetsToLL_01.root -1 1000 2018 MC DY2JetsToLL_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY2JetsToLL_02.root DY2JetsToLL_02.root -1 1000 2018 MC DY2JetsToLL_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY3JetsToLL_00.root DY3JetsToLL_00.root -1 1000 2018 MC DY3JetsToLL_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DY4JetsToLL_00.root DY4JetsToLL_00.root -1 1000 2018 MC DY4JetsToLL_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DYJetsToLL_00.root DYJetsToLL_00.root -1 1000 2018 MC DYJetsToLL_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DYJetsToLL_01.root DYJetsToLL_01.root -1 1000 2018 MC DYJetsToLL_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DYJetsToLL_02.root DYJetsToLL_02.root -1 1000 2018 MC DYJetsToLL_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DYJetsToLL_03.root DYJetsToLL_03.root -1 1000 2018 MC DYJetsToLL_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DYJetsToLL_04.root DYJetsToLL_04.root -1 1000 2018 MC DYJetsToLL_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DYJetsToLL_05.root DYJetsToLL_05.root -1 1000 2018 MC DYJetsToLL_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/DYJetsToLL_06.root DYJetsToLL_06.root -1 1000 2018 MC DYJetsToLL_06 $outDir

./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTTo2L2Nu_powheg_00.root TTTo2L2Nu_powheg_00.root -1 1000 2018 MC TTTo2L2Nu_powheg_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTTo2L2Nu_powheg_01.root TTTo2L2Nu_powheg_01.root -1 1000 2018 MC TTTo2L2Nu_powheg_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTTo2L2Nu_powheg_02.root TTTo2L2Nu_powheg_02.root -1 1000 2018 MC TTTo2L2Nu_powheg_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTTo2L2Nu_powheg_03.root TTTo2L2Nu_powheg_03.root -1 1000 2018 MC TTTo2L2Nu_powheg_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTTo2L2Nu_powheg_04.root TTTo2L2Nu_powheg_04.root -1 1000 2018 MC TTTo2L2Nu_powheg_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTTo2L2Nu_powheg_05.root TTTo2L2Nu_powheg_05.root -1 1000 2018 MC TTTo2L2Nu_powheg_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTTo2L2Nu_powheg_06.root TTTo2L2Nu_powheg_06.root -1 1000 2018 MC TTTo2L2Nu_powheg_06 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_00.root TTToHadronic_powheg_00.root -1 1000 2018 MC TTToHadronic_powheg_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_01.root TTToHadronic_powheg_01.root -1 1000 2018 MC TTToHadronic_powheg_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_02.root TTToHadronic_powheg_02.root -1 1000 2018 MC TTToHadronic_powheg_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_03.root TTToHadronic_powheg_03.root -1 1000 2018 MC TTToHadronic_powheg_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_04.root TTToHadronic_powheg_04.root -1 1000 2018 MC TTToHadronic_powheg_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_05.root TTToHadronic_powheg_05.root -1 1000 2018 MC TTToHadronic_powheg_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_06.root TTToHadronic_powheg_06.root -1 1000 2018 MC TTToHadronic_powheg_06 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_07.root TTToHadronic_powheg_07.root -1 1000 2018 MC TTToHadronic_powheg_07 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_08.root TTToHadronic_powheg_08.root -1 1000 2018 MC TTToHadronic_powheg_08 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_09.root TTToHadronic_powheg_09.root -1 1000 2018 MC TTToHadronic_powheg_09 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_10.root TTToHadronic_powheg_10.root -1 1000 2018 MC TTToHadronic_powheg_10 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_11.root TTToHadronic_powheg_11.root -1 1000 2018 MC TTToHadronic_powheg_11 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_12.root TTToHadronic_powheg_12.root -1 1000 2018 MC TTToHadronic_powheg_12 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_13.root TTToHadronic_powheg_13.root -1 1000 2018 MC TTToHadronic_powheg_13 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToHadronic_powheg_14.root TTToHadronic_powheg_14.root -1 1000 2018 MC TTToHadronic_powheg_14 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToSemiLeptonic_powheg_00.root TTToSemiLeptonic_powheg_00.root -1 1000 2018 MC TTToSemiLeptonic_powheg_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToSemiLeptonic_powheg_01.root TTToSemiLeptonic_powheg_01.root -1 1000 2018 MC TTToSemiLeptonic_powheg_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToSemiLeptonic_powheg_02.root TTToSemiLeptonic_powheg_02.root -1 1000 2018 MC TTToSemiLeptonic_powheg_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToSemiLeptonic_powheg_03.root TTToSemiLeptonic_powheg_03.root -1 1000 2018 MC TTToSemiLeptonic_powheg_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/TTToSemiLeptonic_powheg_04.root TTToSemiLeptonic_powheg_04.root -1 1000 2018 MC TTToSemiLeptonic_powheg_04 $outDir


###########################  DATA #########################
./rootcom etau_analyzer executable_etau  
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_00.root EGamma2018A_00.root -1 1000 2018 DATA EGamma2018A_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_01.root EGamma2018A_01.root -1 1000 2018 DATA EGamma2018A_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_02.root EGamma2018A_02.root -1 1000 2018 DATA EGamma2018A_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_03.root EGamma2018A_03.root -1 1000 2018 DATA EGamma2018A_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_04.root EGamma2018A_04.root -1 1000 2018 DATA EGamma2018A_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_05.root EGamma2018A_05.root -1 1000 2018 DATA EGamma2018A_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_06.root EGamma2018A_06.root -1 1000 2018 DATA EGamma2018A_06 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_07.root EGamma2018A_07.root -1 1000 2018 DATA EGamma2018A_07 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_08.root EGamma2018A_08.root -1 1000 2018 DATA EGamma2018A_08 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_09.root EGamma2018A_09.root -1 1000 2018 DATA EGamma2018A_09 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018A_10.root EGamma2018A_10.root -1 1000 2018 DATA EGamma2018A_10 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018B_00.root EGamma2018B_00.root -1 1000 2018 DATA EGamma2018B_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018B_01.root EGamma2018B_01.root -1 1000 2018 DATA EGamma2018B_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018B_02.root EGamma2018B_02.root -1 1000 2018 DATA EGamma2018B_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018B_03.root EGamma2018B_03.root -1 1000 2018 DATA EGamma2018B_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018B_04.root EGamma2018B_04.root -1 1000 2018 DATA EGamma2018B_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018B_05.root EGamma2018B_05.root -1 1000 2018 DATA EGamma2018B_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018C_00.root EGamma2018C_00.root -1 1000 2018 DATA EGamma2018C_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018C_01.root EGamma2018C_01.root -1 1000 2018 DATA EGamma2018C_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018C_02.root EGamma2018C_02.root -1 1000 2018 DATA EGamma2018C_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018C_03.root EGamma2018C_03.root -1 1000 2018 DATA EGamma2018C_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018C_04.root EGamma2018C_04.root -1 1000 2018 DATA EGamma2018C_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018C_05.root EGamma2018C_05.root -1 1000 2018 DATA EGamma2018C_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018C_06.root EGamma2018C_06.root -1 1000 2018 DATA EGamma2018C_06 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018C_07.root EGamma2018C_07.root -1 1000 2018 DATA EGamma2018C_07 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_00.root EGamma2018D_00.root -1 1000 2018 DATA EGamma2018D_00 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_01.root EGamma2018D_01.root -1 1000 2018 DATA EGamma2018D_01 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_02.root EGamma2018D_02.root -1 1000 2018 DATA EGamma2018D_02 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_03.root EGamma2018D_03.root -1 1000 2018 DATA EGamma2018D_03 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_04.root EGamma2018D_04.root -1 1000 2018 DATA EGamma2018D_04 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_05.root EGamma2018D_05.root -1 1000 2018 DATA EGamma2018D_05 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_06.root EGamma2018D_06.root -1 1000 2018 DATA EGamma2018D_06 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_07.root EGamma2018D_07.root -1 1000 2018 DATA EGamma2018D_07 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_08.root EGamma2018D_08.root -1 1000 2018 DATA EGamma2018D_08 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_09.root EGamma2018D_09.root -1 1000 2018 DATA EGamma2018D_09 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_10.root EGamma2018D_10.root -1 1000 2018 DATA EGamma2018D_10 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_11.root EGamma2018D_11.root -1 1000 2018 DATA EGamma2018D_11 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_12.root EGamma2018D_12.root -1 1000 2018 DATA EGamma2018D_12 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_13.root EGamma2018D_13.root -1 1000 2018 DATA EGamma2018D_13 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_14.root EGamma2018D_14.root -1 1000 2018 DATA EGamma2018D_14 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_15.root EGamma2018D_15.root -1 1000 2018 DATA EGamma2018D_15 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_16.root EGamma2018D_16.root -1 1000 2018 DATA EGamma2018D_16 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_17.root EGamma2018D_17.root -1 1000 2018 DATA EGamma2018D_17 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_18.root EGamma2018D_18.root -1 1000 2018 DATA EGamma2018D_18 $outDir
./MakeCondorFiles.csh executable_etau /hdfs/store/user/jmadhusu/2018_skimmed/etau/EGamma2018D_19.root EGamma2018D_19.root -1 1000 2018 DATA EGamma2018D_19 $outDir
