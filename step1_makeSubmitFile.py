import os
import re


samples = sorted(os.listdir('/hdfs/store/user/jmadhusu/with_boostedtau/2017_skimmed/with_boostedtaus/hadd/'))
outFile = open("submit_condor.sh", "w")
outFile.write("""
./rootcom mutau_analyzer analyze_mutau
outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

""")

for line in samples:
    line = line.strip().replace('.root', '')
    # outFile.write('mkdir -p /hdfs/store/user/jmadhusu/with_boostedtau/2017_skimmed/analyzer/mutau/'+line+' \n')
    if 'SingleMuon' in line:
      line = "./MakeCondorFiles.csh analyze_mutau "+ line + " -1 1000 2017 DATA $outDir \n"
    else:
      line = "./MakeCondorFiles.csh analyze_mutau "+ line + " -1 1000 2017 MC $outDir \n"
  
    outFile.write(line)
    
outFile.close()
print """
created submit_condor.sh
do 
  bash submit_condor.sh

"""
