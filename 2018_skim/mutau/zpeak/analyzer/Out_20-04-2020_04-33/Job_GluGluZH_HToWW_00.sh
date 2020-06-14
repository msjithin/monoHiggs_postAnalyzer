#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc630/cms/cmssw/CMSSW_9_4_13/
eval unset  SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL;
export SRT_CMSSW_VERSION_SCRAMRTDEL="CMSSW_10_2_10";
export SRT_ROOT_INCLUDE_PATH_SCRAMRTDEL="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/src:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/src:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/coral/CORAL_2_3_21-gnimlf9/include/LCG:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/mctester/1.25.0a-gnimlf6/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf5/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/dd4hep/v01-08x-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtDesigner:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/herwigpp/7.1.4-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tauolapp/1.1.5-gnimlf4/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/charybdis/1.003-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/include/SHERPA-MC:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtOpenGL:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtGui:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/thepeg/2.1.4-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia8/230-gnimlf4/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/herwig/6.521-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/rivet/2.5.4-gnimlf6/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/Qt3Support:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lwtnn/2.4-gnimlf4/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/geant4/10.04-gnimlf3/include/Geant4:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/classlib/3.1.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tkonlinesw/4.2.0-1_gcc7-gnimlf6/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/starlight/r193-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/Qt:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtCore:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtXml:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/mcdb/1.0.3-omkpbe2/interface:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libungif/4.1.4-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libtiff/4.0.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libpng/1.6.16-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/frontier_client/2.8.20-omkpbe4/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pcre/8.37-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xrootd/4.8.3-gnimlf/include/xrootd/private:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/vdt/0.4.0-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/valgrind/3.13.0-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/utm/utm_0.7.1-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/toprex/4.23-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tbb/2018_U1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tauola/27.121.5-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sigcpp/2.6.2-omkpbe2/include/sigc++-2.0:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sqlite/3.22.0-omkpbe/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/protobuf/3.5.2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pacparser/1.3.5-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/oracle/12.1.0.2.0/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/meschach/1.2.pCMS1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libuuid/2.22.2-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libhepml/0.2.1-omkpbe2/interface:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/ktjet/1.06-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/jimmy/4.2-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/jemalloc/5.1.0/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/madgraph5amcatnlo/2.6.0-gnimlf8:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/heppdt/3.03.00-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/hector/1.3.4_patch1-gnimlf6/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gsl/2.2.1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libjpeg-turbo/1.3.1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/giflib/4.2.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/vecgeom/v00.05.00-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gdbm/1.10-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/freetype/2.5.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fftw3/3.3.2-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fftjet/1.5.0-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fastjet/3.3.0-omkpbe/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/expat/2.1.0-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/dpm/1.8.0.1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/hepmc/2.06.07-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xerces-c/3.1.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xz/5.2.2-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/dcap/2.47.8-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libxml2/2.9.1-omkpbe2/include/libxml2:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/curl/7.59.0/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cppunit/1.12.1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/clhep/2.4.0.0-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/openssl/1.0.2d-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia6/426-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/photos/215.5-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/zlib-x86_64/1.2.11-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cascade/2.2.04-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/bz2lib/1.0.6-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/python/2.7.14-omkpbe4/include/python2.7:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/ittnotify/16.06.18-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gosamcontrib/2.0-20150803-omkpbe2/include:/usr/local/include:/usr/include";
export SRT_SHERPA_LIBRARY_PATH_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/lib/SHERPA-MC";
export SRT_PYTHON27PATH_SCRAMRTDEL="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/python:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/python:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/coral/CORAL_2_3_21-gnimlf9/slc7_amd64_gcc700/python:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/coral/CORAL_2_3_21-gnimlf9/slc7_amd64_gcc700/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xrootd/4.8.3-gnimlf/lib/python2.7/site-packages:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fastjet/3.3.0-omkpbe/lib/python2.7/site-packages:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw-tool-conf/44.0-gnimlf10/lib/python2.7/site-packages";
export SRT_CMSSW_SEARCH_PATH_SCRAMRTDEL="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/src:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/external/slc7_amd64_gcc700/data:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/src:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/external/slc7_amd64_gcc700/data";
export SRT_PATH_SCRAMRT="/cvmfs/cms.cern.ch/share/overrides/bin:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/bin/slc7_amd64_gcc700:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/external/slc7_amd64_gcc700/bin:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/bin/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/external/slc7_amd64_gcc700/bin:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/llvm/6.0.0-ogkkac/bin:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gcc/7.0.0-omkpbe2/bin";
export SRT_LHAPDF_DATA_PATH_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf2/share/LHAPDF";
export SRT_ROOT_PATH_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/madgraph5amcatnlo/2.6.0-gnimlf8:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gosamcontrib/2.0-20150803-omkpbe2";
export SRT_RIVET_ANALYSIS_PATH_SCRAMRTDEL="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/rivet/2.5.4-gnimlf6/lib";
export SRT_THEPEGPATH_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/thepeg/2.1.4-gnimlf/share/ThePEG";
export SRT_CASCADE_PDFPATH_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cascade/2.2.04-gnimlf2/share";
export SRT_VINCIADATA_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/vincia/2.2.01-gnimlf4/share/Vincia/xmldoc";
export SRT_CMSSW_BASE_SCRAMRTDEL="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10";
export SRT_HERWIGPATH_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/herwigpp/7.1.4-gnimlf/share/Herwig";
export SRT_EVTGENDATA_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/evtgen/1.6.0-gnimlf4/share";
export SRT_LOCALRT_SCRAMRTDEL="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10";
export SRT_SHERPA_SHARE_PATH_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/share/SHERPA-MC";
export SRT_SHERPA_INCLUDE_PATH_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/include/SHERPA-MC";
export SCRAMRT_SET="CMSSW:CMSSW_10_2_10:slc7_amd64_gcc700:V2_2_9_pre06:SRT_";
export SRT_LD_LIBRARY_PATH_SCRAMRTDEL="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/biglib/slc7_amd64_gcc700:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/external/slc7_amd64_gcc700/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/biglib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/external/slc7_amd64_gcc700/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/llvm/6.0.0-ogkkac/lib64:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gcc/7.0.0-omkpbe2/lib64:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gcc/7.0.0-omkpbe2/lib";
export SRT_PYTHON3PATH_SCRAMRTDEL="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/python:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/python:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw-tool-conf/44.0-gnimlf10/lib/python3.6/site-packages";
export SRT_PYTHIA8DATA_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia8/230-gnimlf4/share/Pythia8/xmldoc";
export SRT_CMSSW_GIT_HASH_SCRAMRTDEL="CMSSW_10_2_10";
export SRT_CMSSW_FWLITE_INCLUDE_PATH_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/hepmc/2.06.07-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/clhep/2.4.0.0-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0-gnimlf/include";
export SRT_CMSSW_RELEASE_BASE_SCRAMRT="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10";
export SRT_CMS_OPENLOOPS_PREFIX_SCRAMRTDEL="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/openloops/2.0.b";
export CMSSW_BASE="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10";
export CMSSW_VERSION="CMSSW_10_2_10";
export CMSSW_GIT_HASH="CMSSW_10_2_10";
export LOCALRT="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10";
export CMSSW_RELEASE_BASE="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10";
export EVTGENDATA="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/evtgen/1.6.0-gnimlf4/share";
export VINCIADATA="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/vincia/2.2.01-gnimlf4/share/Vincia/xmldoc";
export HERWIGPATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/herwigpp/7.1.4-gnimlf/share/Herwig";
export THEPEGPATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/thepeg/2.1.4-gnimlf/share/ThePEG";
export PYTHIA8DATA="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia8/230-gnimlf4/share/Pythia8/xmldoc";
export LHAPDF_DATA_PATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf2/share/LHAPDF";
export CASCADE_PDFPATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cascade/2.2.04-gnimlf2/share";
export CMSSW_FWLITE_INCLUDE_PATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/hepmc/2.06.07-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/clhep/2.4.0.0-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0-gnimlf/include";
export SHERPA_LIBRARY_PATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/lib/SHERPA-MC";
export SHERPA_INCLUDE_PATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/include/SHERPA-MC";
export PATH="/cvmfs/cms.cern.ch/share/overrides/bin:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/bin/slc7_amd64_gcc700:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/external/slc7_amd64_gcc700/bin:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/bin/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/external/slc7_amd64_gcc700/bin:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/llvm/6.0.0-ogkkac/bin:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gcc/7.0.0-omkpbe2/bin:/cvmfs/cms.cern.ch/common:/cvmfs/cms.cern.ch/common:/usr/local/bin/smb:/usr/local/bin/raid:/usr/lib64/qt-3.3/bin:/usr/local/etc:/usr/local/bin/hadoop:/usr/local/bin/firstboot:/usr/local/bin/cephdir:/usr/local/bin/afs:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/afs/hep.wisc.edu/home/ms/bin";
export CMSSW_SEARCH_PATH="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/src:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/external/slc7_amd64_gcc700/data:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/src:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/external/slc7_amd64_gcc700/data";
export ROOT_PATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/madgraph5amcatnlo/2.6.0-gnimlf8:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gosamcontrib/2.0-20150803-omkpbe2";
export RIVET_ANALYSIS_PATH="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/rivet/2.5.4-gnimlf6/lib";
export LD_LIBRARY_PATH="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/biglib/slc7_amd64_gcc700:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/external/slc7_amd64_gcc700/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/biglib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/external/slc7_amd64_gcc700/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/llvm/6.0.0-ogkkac/lib64:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gcc/7.0.0-omkpbe2/lib64:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gcc/7.0.0-omkpbe2/lib";
export CMS_OPENLOOPS_PREFIX="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/openloops/2.0.b";
export ROOT_INCLUDE_PATH="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/src:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/src:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/coral/CORAL_2_3_21-gnimlf9/include/LCG:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/mctester/1.25.0a-gnimlf6/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf5/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/dd4hep/v01-08x-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtDesigner:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/herwigpp/7.1.4-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tauolapp/1.1.5-gnimlf4/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/charybdis/1.003-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/include/SHERPA-MC:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtOpenGL:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtGui:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/thepeg/2.1.4-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia8/230-gnimlf4/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/herwig/6.521-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/rivet/2.5.4-gnimlf6/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/Qt3Support:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lwtnn/2.4-gnimlf4/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/geant4/10.04-gnimlf3/include/Geant4:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/classlib/3.1.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tkonlinesw/4.2.0-1_gcc7-gnimlf6/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/starlight/r193-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/Qt:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtCore:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtXml:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/mcdb/1.0.3-omkpbe2/interface:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libungif/4.1.4-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libtiff/4.0.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libpng/1.6.16-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/frontier_client/2.8.20-omkpbe4/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pcre/8.37-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xrootd/4.8.3-gnimlf/include/xrootd/private:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/vdt/0.4.0-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/valgrind/3.13.0-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/utm/utm_0.7.1-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/toprex/4.23-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tbb/2018_U1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tauola/27.121.5-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sigcpp/2.6.2-omkpbe2/include/sigc++-2.0:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sqlite/3.22.0-omkpbe/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/protobuf/3.5.2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pacparser/1.3.5-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/oracle/12.1.0.2.0/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/meschach/1.2.pCMS1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libuuid/2.22.2-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libhepml/0.2.1-omkpbe2/interface:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/ktjet/1.06-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/jimmy/4.2-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/jemalloc/5.1.0/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/madgraph5amcatnlo/2.6.0-gnimlf8:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/heppdt/3.03.00-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/hector/1.3.4_patch1-gnimlf6/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gsl/2.2.1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libjpeg-turbo/1.3.1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/giflib/4.2.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/vecgeom/v00.05.00-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gdbm/1.10-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/freetype/2.5.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fftw3/3.3.2-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fftjet/1.5.0-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fastjet/3.3.0-omkpbe/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/expat/2.1.0-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/dpm/1.8.0.1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/hepmc/2.06.07-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xerces-c/3.1.3-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xz/5.2.2-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/dcap/2.47.8-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libxml2/2.9.1-omkpbe2/include/libxml2:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/curl/7.59.0/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cppunit/1.12.1-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/clhep/2.4.0.0-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/openssl/1.0.2d-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia6/426-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/photos/215.5-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/zlib-x86_64/1.2.11-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cascade/2.2.04-gnimlf2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/bz2lib/1.0.6-omkpbe2/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/python/2.7.14-omkpbe4/include/python2.7:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/ittnotify/16.06.18-gnimlf/include:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gosamcontrib/2.0-20150803-omkpbe2/include:/usr/local/include:/usr/include";
export PYTHON3PATH="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/python:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/python:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw-tool-conf/44.0-gnimlf10/lib/python3.6/site-packages";
export SHERPA_SHARE_PATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.5-gnimlf/share/SHERPA-MC";
export PYTHON27PATH="/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/python:/afs/hep.wisc.edu/user/ms/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/python:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_10/lib/slc7_amd64_gcc700:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/coral/CORAL_2_3_21-gnimlf9/slc7_amd64_gcc700/python:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/coral/CORAL_2_3_21-gnimlf9/slc7_amd64_gcc700/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xrootd/4.8.3-gnimlf/lib/python2.7/site-packages:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fastjet/3.3.0-omkpbe/lib/python2.7/site-packages:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw-tool-conf/44.0-gnimlf10/lib/python2.7/site-packages";
cd ${_CONDOR_SCRATCH_DIR}
./analyze_mutau /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/GluGluZH_HToWW_00.root GluGluZH_HToWW_00.root -1 1000 2018_test MC GluGluZH_HToWW_00

