
./rootcom skimm_tt_2018 analyze_tautau_skim 
outDir="Out_DATA_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauA/200609_122655/0000/ TauA_00.root -1 1000 2018 DATA TauA_00 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauA/200609_122655/0001/ TauA_01.root -1 1000 2018 DATA TauA_01 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauA/200609_122655/0002/ TauA_02.root -1 1000 2018 DATA TauA_02 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauB/200609_122722/0000/ TauB_00.root -1 1000 2018 DATA TauB_00 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauB/200609_122722/0001/ TauB_01.root -1 1000 2018 DATA TauB_01 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauB/200609_122722/0002/ TauB_02.root -1 1000 2018 DATA TauB_02 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauC/200609_122756/0000/ TauC_00.root -1 1000 2018 DATA TauC_00 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauC/200609_122756/0001/ TauC_01.root -1 1000 2018 DATA TauC_01 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauC/200609_122756/0002/ TauC_02.root -1 1000 2018 DATA TauC_02 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0000/ TauD_PromptReco_00.root -1 1000 2018 DATA TauD_PromptReco_00 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0001/ TauD_PromptReco_01.root -1 1000 2018 DATA TauD_PromptReco_01 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0002/ TauD_PromptReco_02.root -1 1000 2018 DATA TauD_PromptReco_02 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0003/ TauD_PromptReco_03.root -1 1000 2018 DATA TauD_PromptReco_03 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0004/ TauD_PromptReco_04.root -1 1000 2018 DATA TauD_PromptReco_04 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0005/ TauD_PromptReco_05.root -1 1000 2018 DATA TauD_PromptReco_05 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0006/ TauD_PromptReco_06.root -1 1000 2018 DATA TauD_PromptReco_06 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0007/ TauD_PromptReco_07.root -1 1000 2018 DATA TauD_PromptReco_07 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0008/ TauD_PromptReco_08.root -1 1000 2018 DATA TauD_PromptReco_08 $outDir 
./MakeCondorFiles.csh analyze_tautau_skim /hdfs/store/user/jmadhusu/data2018_09Jun2020/Tau/crab_job_TauD_PromptReco/200609_122856/0009/ TauD_PromptReco_09.root -1 1000 2018 DATA TauD_PromptReco_09 $outDir 
