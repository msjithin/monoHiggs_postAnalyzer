
./rootcom skimm_mt_2017 analyze_mutau_skim 
outDir="Out_MC_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_job_TTToSemiLeptonic_TuneCP5/200609_002140/0000/ TTToSemiLeptonic_TuneCP5_00.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_job_TTToSemiLeptonic_TuneCP5/200609_002140/0001/ TTToSemiLeptonic_TuneCP5_01.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_job_TTToSemiLeptonic_TuneCP5/200609_002140/0002/ TTToSemiLeptonic_TuneCP5_02.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_job_TTToSemiLeptonic_TuneCP5/200609_002140/0003/ TTToSemiLeptonic_TuneCP5_03.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_03 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_job_TTToSemiLeptonic_TuneCP5/200609_002140/0004/ TTToSemiLeptonic_TuneCP5_04.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_04 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_job_TTToSemiLeptonic_TuneCP5/200609_002140/0005/ TTToSemiLeptonic_TuneCP5_05.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_05 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_job_TTToSemiLeptonic_TuneCP5/200609_002140/0006/ TTToSemiLeptonic_TuneCP5_06.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_06 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_job_TTToSemiLeptonic_TuneCP5/200609_002140/0007/ TTToSemiLeptonic_TuneCP5_07.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_07 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_09Jun2020/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_job_TTToSemiLeptonic_TuneCP5/200609_002140/0008/ TTToSemiLeptonic_TuneCP5_08.root -1 1000 2017 MC TTToSemiLeptonic_TuneCP5_08 $outDir 
