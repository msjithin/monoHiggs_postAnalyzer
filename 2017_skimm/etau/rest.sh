
./rootcom skimm_et_2017 analyze_etau_skim 
outDir="Out_MC_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

./MakeCondorFiles.csh analyze_etau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_job_TTTo2L2Nu_TuneCP5/200609_002055/0000/ TTTo2L2Nu_TuneCP5_00.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_00 $outDir 
./MakeCondorFiles.csh analyze_etau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_job_TTTo2L2Nu_TuneCP5/200609_002055/0001/ TTTo2L2Nu_TuneCP5_01.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_01 $outDir 
./MakeCondorFiles.csh analyze_etau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_job_TTTo2L2Nu_TuneCP5/200609_002055/0002/ TTTo2L2Nu_TuneCP5_02.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_02 $outDir 
./MakeCondorFiles.csh analyze_etau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_job_TTTo2L2Nu_TuneCP5/200609_002055/0003/ TTTo2L2Nu_TuneCP5_03.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_03 $outDir 
./MakeCondorFiles.csh analyze_etau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_job_TTTo2L2Nu_TuneCP5/200609_002055/0004/ TTTo2L2Nu_TuneCP5_04.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_04 $outDir 
./MakeCondorFiles.csh analyze_etau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_job_TTTo2L2Nu_TuneCP5/200609_002055/0005/ TTTo2L2Nu_TuneCP5_05.root -1 1000 2017 MC TTTo2L2Nu_TuneCP5_05 $outDir 
