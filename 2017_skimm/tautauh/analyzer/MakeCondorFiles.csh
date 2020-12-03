#!/bin/sh

currentDir=${PWD}
outDir=${9}
#export CMSSW_RELEASE_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_1/
export CMSSW_RELEASE_BASE=/cvmfs/cms.cern.ch/slc7_amd64_gcc630/cms/cmssw/CMSSW_9_4_13/
cat> $outDir/Job_${8}.sh<<EOF
#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd $CMSSW_RELEASE_BASE
eval `scramv1 runtime -sh`
cd \${_CONDOR_SCRATCH_DIR}
./${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8}

EOF

chmod 775 $outDir/Job_${8}.sh
cat > $outDir/condor_${8}<<EOF
x509userproxy = /tmp/x509up_u4548
universe = vanilla
Executable = $outDir/Job_${8}.sh
Notification         = never
WhenToTransferOutput = On_Exit
ShouldTransferFiles  = yes
Requirements = (TARGET.UidDomain == "hep.wisc.edu" && TARGET.HAS_CMS_HDFS && OpSysAndVer == "CENTOS7" && TARGET.Arch == "X86_64" && (MY.RequiresSharedFS=!=true || TARGET.HasAFS_OSG) && (TARGET.OSG_major =!= undefined || TARGET.IS_GLIDEIN=?=true) && IsSlowSlot=!=true)
on_exit_remove       = (ExitBySignal == FALSE && (ExitCode == 0 || ExitCode == 42 || NumJobStarts>3))
+IsFastQueueJob      = True
getenv               = True
request_memory       = 1992
request_disk         = 2048000

#OutputDestination = ${outdir}
#Initialdir = Out_${8}         
Transfer_Input_Files = ${currentDir}/${1} , /nfs_scratch/jmadhusu/CMSSW_10_2_10/src/sf_files

output               = $outDir/\$(Cluster)_\$(Process)_${8}.out
error                = $outDir/\$(Cluster)_\$(Process)_${8}.err
Log                  = $outDir/\$(Cluster)_\$(Process)_${8}.log
Queue
EOF

condor_submit $outDir/condor_${8}
