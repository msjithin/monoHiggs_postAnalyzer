#!/bin/sh

currentDir=${PWD}
outDir=${7}
#export CMSSW_RELEASE_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_1/
export CMSSW_RELEASE_BASE=/cvmfs/cms.cern.ch/slc7_amd64_gcc630/cms/cmssw/CMSSW_9_4_13/
cat> $outDir/Job_${2}.sh<<EOF
#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd $CMSSW_RELEASE_BASE
eval `scramv1 runtime -sh`
cd \${_CONDOR_SCRATCH_DIR}

start_time=\$(date +%s)

count=0

./${1} root://cmsxrootd.hep.wisc.edu:1094//store/user/jmadhusu/with_boostedtau/2017_skimmed/with_boostedtaus/hadd/${2}.root ${2}.root ${3} ${4} ${5} ${6} ${2}

end_time=\$(date +%s)

duration=\$((\$end_time-\$start_time))
echo "Total time \$((\$duration / 3600)) hours , \$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed."

EOF

chmod 775 $outDir/Job_${2}.sh
cat > $outDir/condor_${2}<<EOF
x509userproxy = /tmp/x509up_u4548
universe = vanilla
Executable = $outDir/Job_${2}.sh
Notification         = never
WhenToTransferOutput = On_Exit
ShouldTransferFiles  = yes
Requirements = ( OpSysAndVer == "CENTOS7" && TARGET.Arch == "X86_64" && (MY.RequiresSharedFS=!=true || TARGET.HasAFS_OSG) && (TARGET.HAS_OSG_WN_CLIENT =?= TRUE || TARGET.IS_GLIDEIN=?=true) && IsSlowSlot=!=true)
on_exit_remove       = (ExitBySignal == FALSE && (ExitCode == 0 || ExitCode == 42 || NumJobStarts>3))
+IsFastQueueJob      = True
getenv               = True
request_memory       = 8G
request_disk         = 8G

#OutputDestination = ${outdir}
#Initialdir = Out_${7}         
Transfer_Input_Files = ${currentDir}/${1} , /nfs_scratch/jmadhusu/CMSSW_10_2_10/src/sf_files 

output               = $outDir/\$(Cluster)_\$(Process)_${2}.out
error                = $outDir/\$(Cluster)_\$(Process)_${2}.err
Log                  = $outDir/\$(Cluster)_\$(Process)_${2}.log
Queue 
EOF

condor_submit $outDir/condor_${2}
