

#./rootcom Analyzer_mutau_data analyze_mutau_data
#./analyze_mutau_data /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraB/200426_094119/0000/ SingleMuon_test.root 10000 10000 2017 DATA SingleMuon_10

./rootcom calc_FR_2017 analyze_mutau_updated
./analyze_mutau_updated /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/crab_job_SingleMuon_EraB/200426_094119/0000/ SingleMuon_test.root 100000 10000 2017 DATA

#./rootcom Analyzer_etau_data analyze_etau_data
#./analyze_etau_data /hdfs/store/user/gomber/MonoHiggs_2017_Data/monoHiggs_data_2017_31March2018/SingleElectron/crab_sEle_dataset1/181008_092222/0000/ SingleElectron_10_test.root 1000000 100000 SingleElectron_10
