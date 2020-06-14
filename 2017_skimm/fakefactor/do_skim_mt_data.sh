
./rootcom calc_FR_2017 analyze_mutau_skim 
outDir="Out_DATA_$(date +"%d-%m-%Y_%H-%M")" 
outName="$(date +"%d-%m-%Y_%H-%M")_"
mkdir $outDir 

cp calc_FR_2017.C $outDir
cp calc_FR_2017.h $outDir


./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraB/200426_094119/0000/ ${outName}SingleMuon_EraB_00.root -1 1000 2017 DATA SingleMuon_EraB_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraB/200426_094119/0001/ ${outName}SingleMuon_EraB_01.root -1 1000 2017 DATA SingleMuon_EraB_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraB/200426_094119/0002/ ${outName}SingleMuon_EraB_02.root -1 1000 2017 DATA SingleMuon_EraB_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraB/200426_094119/0003/ ${outName}SingleMuon_EraB_03.root -1 1000 2017 DATA SingleMuon_EraB_03 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraB/200426_094119/0004/ ${outName}SingleMuon_EraB_04.root -1 1000 2017 DATA SingleMuon_EraB_04 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraB/200426_094119/0005/ ${outName}SingleMuon_EraB_05.root -1 1000 2017 DATA SingleMuon_EraB_05 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraC/200426_094141/0000/ ${outName}SingleMuon_EraC_00.root -1 1000 2017 DATA SingleMuon_EraC_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraC/200426_094141/0001/ ${outName}SingleMuon_EraC_01.root -1 1000 2017 DATA SingleMuon_EraC_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraC/200426_094141/0002/ ${outName}SingleMuon_EraC_02.root -1 1000 2017 DATA SingleMuon_EraC_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraC/200426_094141/0003/ ${outName}SingleMuon_EraC_03.root -1 1000 2017 DATA SingleMuon_EraC_03 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraC/200426_094141/0004/ ${outName}SingleMuon_EraC_04.root -1 1000 2017 DATA SingleMuon_EraC_04 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraD/200426_094207/0000/ ${outName}SingleMuon_EraD_00.root -1 1000 2017 DATA SingleMuon_EraD_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraD/200426_094207/0001/ ${outName}SingleMuon_EraD_01.root -1 1000 2017 DATA SingleMuon_EraD_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraD/200426_094207/0002/ ${outName}SingleMuon_EraD_02.root -1 1000 2017 DATA SingleMuon_EraD_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraE/200426_094230/0000/ ${outName}SingleMuon_EraE_00.root -1 1000 2017 DATA SingleMuon_EraE_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraE/200426_094230/0001/ ${outName}SingleMuon_EraE_01.root -1 1000 2017 DATA SingleMuon_EraE_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraE/200426_094230/0002/ ${outName}SingleMuon_EraE_02.root -1 1000 2017 DATA SingleMuon_EraE_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraE/200426_094230/0003/ ${outName}SingleMuon_EraE_03.root -1 1000 2017 DATA SingleMuon_EraE_03 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraE/200426_094230/0004/ ${outName}SingleMuon_EraE_04.root -1 1000 2017 DATA SingleMuon_EraE_04 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraE/200426_094230/0005/ ${outName}SingleMuon_EraE_05.root -1 1000 2017 DATA SingleMuon_EraE_05 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraF/200426_094254/0000/ ${outName}SingleMuon_EraF_00.root -1 1000 2017 DATA SingleMuon_EraF_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraF/200426_094254/0001/ ${outName}SingleMuon_EraF_01.root -1 1000 2017 DATA SingleMuon_EraF_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraF/200426_094254/0002/ ${outName}SingleMuon_EraF_02.root -1 1000 2017 DATA SingleMuon_EraF_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraF/200426_094254/0003/ ${outName}SingleMuon_EraF_03.root -1 1000 2017 DATA SingleMuon_EraF_03 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraF/200426_094254/0004/ ${outName}SingleMuon_EraF_04.root -1 1000 2017 DATA SingleMuon_EraF_04 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraF/200426_094254/0005/ ${outName}SingleMuon_EraF_05.root -1 1000 2017 DATA SingleMuon_EraF_05 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraF/200426_094254/0006/ ${outName}SingleMuon_EraF_06.root -1 1000 2017 DATA SingleMuon_EraF_06 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraF/200426_094254/0007/ ${outName}SingleMuon_EraF_07.root -1 1000 2017 DATA SingleMuon_EraF_07 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraF/200426_094254/0008/ ${outName}SingleMuon_EraF_08.root -1 1000 2017 DATA SingleMuon_EraF_08 $outDir 
