
./rootcom skimm_et_2017 analyze_etau_skim 
outDir="Out_DATA_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

./MakeCondorFiles.csh analyze_etau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleElectron/crab_job_SingleElectron_EraB/200426_093924/0000/ SingleElectron_EraB_00.root -1 1000 2017 DATA SingleElectron_EraB_00 $outDir 
./MakeCondorFiles.csh analyze_etau_skim /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleElectron/crab_job_SingleElectron_EraB/200426_093924/0001/ SingleElectron_EraB_01.root -1 1000 2017 DATA SingleElectron_EraB_01 $outDir 

