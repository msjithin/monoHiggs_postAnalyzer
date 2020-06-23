
./rootcom skimm_mt_2018 analyze_mutau_skim 
outDir="Out_DATA_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonA/200609_122536/0000/ SingleMuonA_00.root -1 1000 2018 DATA SingleMuonA_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonA/200609_122536/0001/ SingleMuonA_01.root -1 1000 2018 DATA SingleMuonA_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonA/200609_122536/0002/ SingleMuonA_02.root -1 1000 2018 DATA SingleMuonA_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonA/200609_122536/0003/ SingleMuonA_03.root -1 1000 2018 DATA SingleMuonA_03 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonA/200609_122536/0004/ SingleMuonA_04.root -1 1000 2018 DATA SingleMuonA_04 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonA/200609_122536/0005/ SingleMuonA_05.root -1 1000 2018 DATA SingleMuonA_05 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonA/200609_122536/0006/ SingleMuonA_06.root -1 1000 2018 DATA SingleMuonA_06 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonB/200609_122603/0000/ SingleMuonB_00.root -1 1000 2018 DATA SingleMuonB_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonB/200609_122603/0001/ SingleMuonB_01.root -1 1000 2018 DATA SingleMuonB_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonB/200609_122603/0002/ SingleMuonB_02.root -1 1000 2018 DATA SingleMuonB_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonB/200609_122603/0003/ SingleMuonB_03.root -1 1000 2018 DATA SingleMuonB_03 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonB/200609_122603/0004/ SingleMuonB_04.root -1 1000 2018 DATA SingleMuonB_04 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonC/200609_122629/0000/ SingleMuonC_00.root -1 1000 2018 DATA SingleMuonC_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonC/200609_122629/0001/ SingleMuonC_01.root -1 1000 2018 DATA SingleMuonC_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonC/200609_122629/0002/ SingleMuonC_02.root -1 1000 2018 DATA SingleMuonC_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonC/200609_122629/0003/ SingleMuonC_03.root -1 1000 2018 DATA SingleMuonC_03 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0000/ SingleMuonD_PromptReco_00.root -1 1000 2018 DATA SingleMuonD_PromptReco_00 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0001/ SingleMuonD_PromptReco_01.root -1 1000 2018 DATA SingleMuonD_PromptReco_01 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0002/ SingleMuonD_PromptReco_02.root -1 1000 2018 DATA SingleMuonD_PromptReco_02 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0003/ SingleMuonD_PromptReco_03.root -1 1000 2018 DATA SingleMuonD_PromptReco_03 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0004/ SingleMuonD_PromptReco_04.root -1 1000 2018 DATA SingleMuonD_PromptReco_04 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0005/ SingleMuonD_PromptReco_05.root -1 1000 2018 DATA SingleMuonD_PromptReco_05 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0006/ SingleMuonD_PromptReco_06.root -1 1000 2018 DATA SingleMuonD_PromptReco_06 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0007/ SingleMuonD_PromptReco_07.root -1 1000 2018 DATA SingleMuonD_PromptReco_07 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0008/ SingleMuonD_PromptReco_08.root -1 1000 2018 DATA SingleMuonD_PromptReco_08 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0009/ SingleMuonD_PromptReco_09.root -1 1000 2018 DATA SingleMuonD_PromptReco_09 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0010/ SingleMuonD_PromptReco_10.root -1 1000 2018 DATA SingleMuonD_PromptReco_10 $outDir 
./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/SingleMuon/crab_job_SingleMuonD_PromptReco/200609_122829/0011/ SingleMuonD_PromptReco_11.root -1 1000 2018 DATA SingleMuonD_PromptReco_11 $outDir 
