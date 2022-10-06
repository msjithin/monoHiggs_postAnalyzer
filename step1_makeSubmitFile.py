import os
import re

os.popen('ls -d /hdfs/store/user/jmadhusu/with_boostedtau/2017_skimmed/with_boostedtaus/mutau/* > mc_dir_list')
inputFile=open("mc_dir_list", "r")
outFile = open("submit_condor.sh", "w")
outFile.write("""
./rootcom mutau_analyzer analyze_mutau
outDir="Out_MC_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

""")

for line in inputFile:
    line = line.strip()
    line = line.replace( "/hdfs/store/","./MakeCondorFiles.csh analyze_mutau /store/")
    sample_name = line.split('/')[-1]
    print sample_name
    if 'SingleMuon' in line:
      line = line + ' ' + sample_name + ' -1 1000 2017 DATA ' + sample_name + ' $outDir \n'
    else:
      line = line + ' ' + sample_name + ' -1 1000 2017 MC ' + sample_name + ' $outDir \n'
    outFile.write(line)

inputFile.close()
outFile.close()
print """
created submit_condor.sh
do 
  bash submit_condor.sh

"""
