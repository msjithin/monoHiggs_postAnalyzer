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

start_time=\$(date +%s)

xrdfs cmsxrootd.hep.wisc.edu ls ${2} >> filelist

count=0
while IFS= read -r line
do
  echo \$line
  count=\$((count+1))
  ./${1} root://cmsxrootd.hep.wisc.edu:1094/\$line ${8}_\${count}.root ${4} ${5} ${6} ${7} ${8} &

  #xrdcp ${8}_\${count}.root root://cmsxrootd.hep.wisc.edu:1094//store/user/jmadhusu/with_boostedtau/2017_skimmed/with_boostedtaus/tautau/${8}/
  #rm ${8}_\${count}.root 
  if (( \$count % 5 == 0 ))
  then
    sleep 30
  fi

done < "filelist"

# sleep 60

# hadd ${8}.root ${8}_*.root
# rm  ${8}_*.root
rm filelist

end_time=\$(date +%s)

duration=\$((\$end_time-\$start_time))
echo "Total time \$((\$duration / 3600)) hours , \$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed."

EOF

chmod 775 $outDir/Job_${8}.sh
cat > $outDir/condor_${8}<<EOF
x509userproxy = /tmp/x509up_u4548
universe = vanilla
Executable = $outDir/Job_${8}.sh
Notification         = never
WhenToTransferOutput = On_Exit
ShouldTransferFiles  = yes
Requirements = ( OpSysAndVer == "CENTOS7" && TARGET.Arch == "X86_64" && (MY.RequiresSharedFS=!=true || TARGET.HasAFS_OSG) && (TARGET.HAS_OSG_WN_CLIENT =?= TRUE || TARGET.IS_GLIDEIN=?=true) && IsSlowSlot=!=true)
on_exit_remove       = (ExitBySignal == FALSE && (ExitCode == 0 || ExitCode == 42 || NumJobStarts>3))
+IsFastQueueJob      = True
getenv               = True
request_memory       = 10G
request_disk         = 10G

#OutputDestination = ${outdir}
#Initialdir = Out_${8}         
Transfer_Input_Files = ${currentDir}/${1} , /nfs_scratch/jmadhusu/CMSSW_10_2_10/src/sf_files

output               = $outDir/\$(Cluster)_\$(Process)_${8}.out
error                = $outDir/\$(Cluster)_\$(Process)_${8}.err
Log                  = $outDir/\$(Cluster)_\$(Process)_${8}.log
Queue
EOF

condor_submit $outDir/condor_${8}
