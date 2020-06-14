
./rootcom Analyzer_mutau_data analyze_mutau_data
outDir="Out_DATA_$(date +"%d-%m-%Y_%H-%M")" 
outfile="SingleMuon_$(date +"%d-%m-%Y_%H-%M").root"
mkdir $outDir 
cp Analyzer_mutau_data* $outDir


#./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset1/181008_091856/0000/ "10_"$outfile 1000 1000 2017 DATA SingleMuon_10 $outDir 

./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset1/181008_091856/0000/ "10_"$outfile -1 1000000 2017 DATA SingleMuon_10 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset1/181008_091856/0001/ "11_"$outfile -1 1000000 2017 DATA SingleMuon_11 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset2/181008_091921/0000/ "20_"$outfile -1 1000000 2017 DATA SingleMuon_20 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset2/181008_091921/0001/ "21_"$outfile -1 1000000 2017 DATA SingleMuon_21 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset2/181008_091921/0002/ "22_"$outfile -1 1000000 2017 DATA SingleMuon_22 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset2/181008_091921/0003/ "23_"$outfile -1 1000000 2017 DATA SingleMuon_23 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset3/181008_091946/0000/ "30_"$outfile -1 1000000 2017 DATA SingleMuon_30 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset3/181008_091946/0001/ "31_"$outfile -1 1000000 2017 DATA SingleMuon_31 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset4/181008_092008/0000/ "40_"$outfile -1 1000000 2017 DATA SingleMuon_40 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset4/181008_092008/0001/ "41_"$outfile -1 1000000 2017 DATA SingleMuon_41 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset4/181008_092008/0002/ "42_"$outfile -1 1000000 2017 DATA SingleMuon_42 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset5/181008_092031/0000/ "50_"$outfile -1 1000000 2017 DATA SingleMuon_50 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset5/181008_092031/0001/ "51_"$outfile -1 1000000 2017 DATA SingleMuon_51 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset5/181008_092031/0002/ "52_"$outfile -1 1000000 2017 DATA SingleMuon_52 $outDir
./MakeCondorFiles.csh analyze_mutau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleMuon/crab_sMu_dataset5/181008_092031/0003/ "53_"$outfile -1 1000000 2017 DATA SingleMuon_53 $outDir
