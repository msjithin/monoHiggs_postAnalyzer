#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc630/cms/cmssw/CMSSW_9_4_13/
eval unset  SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL;
export SRT_CMSSW_RELEASE_BASE_SCRAMRT="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10";
export CMSSW_RELEASE_BASE="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10";
export PATH="/cvmfs/cms.cern.ch/share/overrides/bin:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/bin/slc7_amd64_gcc700:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/external/slc7_amd64_gcc700/bin:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/bin/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/external/slc7_amd64_gcc700/bin:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/llvm/6.0.0-ogkkac/bin:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gcc/7.0.0-omkpbe2/bin:/cvmfs/cms.cern.ch/common:/cvmfs/cms.cern.ch/common:/usr/local/bin/smb:/usr/local/bin/raid:/usr/lib64/qt-3.3/bin:/usr/local/etc:/usr/local/bin/hadoop:/usr/local/bin/firstboot:/usr/local/bin/cephdir:/usr/local/bin/afs:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/afs/hep.wisc.edu/home/ms/bin";
cd ${_CONDOR_SCRATCH_DIR}
./analyze_mutau /hdfs/store/user/jmadhusu/2017_skimmed/mutau/DYJetsToLL_M-50_TuneCP5_v1_02.root DYJetsToLL_M-50_TuneCP5_v1_02.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_02

