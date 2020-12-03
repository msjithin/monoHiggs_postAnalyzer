#!/bin/bash
set -e  # exit when any command fails

nEvents=-1
sample='dy'
plottingOn=0
while getopts n:s:p option
do
    case "${option}"
	in
	n) nEvents=${OPTARG};;
	p) plottingOn=1 ;;
	s) sample=${OPTARG};;
esac
done

#dir :  /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/
#dir : /hdfs/store/user/jmadhusu/fromCecile/skimmedFSA/tautau/
./rootcom smhtt_2017 executable_smhtt_tt
./executable_smhtt_tt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/DY.root theirs_rootfile/DY.root $nEvents 1000 2017 MC DY_v1
#./executable_smhtt_tt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/DY1.root theirs_rootfile/DY1.root $nEvents 1000 2017 MC DY1_v1
#./executable_smhtt_tt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/DY2.root theirs_rootfile/DY2.root $nEvents 1000 2017 MC DY2_v1
#./executable_smhtt_tt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/DY3.root theirs_rootfile/DY3.root $nEvents 1000 2017 MC DY3_v1
#./executable_smhtt_tt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/DY4.root theirs_rootfile/DY4.root $nEvents 1000 2017 MC DY4_v1
./executable_smhtt_tt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/TTTo2L2Nu.root theirs_rootfile/TTTo2L2Nu.root $nEvents 1000 2017 MC TTTo2L2Nu
./executable_smhtt_tt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/TTToHadronic.root theirs_rootfile/TTToHadronic.root $nEvents 1000 2017 MC TTToHadronic
#./executable_smhtt_tt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/TTToSemiLeptonic.root theirs_rootfile/TTToSemiLeptonic.root $nEvents 1000 2017 MC TTToSemiLeptonic



./executable_smhtt_tt /hdfs/store/user/jmadhusu/fromCecile/svfitted_20nov/tautau2017/Data.root theirs_rootfile/Data.root $nEvents 1000 2017 DATA Data
