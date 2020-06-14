./rootcom skimm_mt_2018 analyze_mutau_skim  

outDir="Out_MC_$(date +"%d-%m-%Y_%H-%M")"
mkdir $outDir


./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/user/jmadhusu/MC2018_Autumn18_monoHiggs_28Mar2020/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_DY1JetsToLL/200328_192020/0000/ DY1JetsToLL_00.root 10000 1000 2018_test MC DY1JetsToLL_00 $outDir
