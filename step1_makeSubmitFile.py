import os
import re


samples = sorted(os.listdir('/hdfs/store/user/jmadhusu/with_boostedtau/2017_skimmed/with_boostedtaus/hadd/'))

# crate data and mc list in samples.txt
outFile = open("samples.txt", "w")

for line in samples:
    line = line.strip().replace('.root', '')
    # outFile.write('mkdir -p /hdfs/store/user/jmadhusu/with_boostedtau/2017_skimmed/analyzer/mutau/'+line+' \n')
    if 'SingleMuon' in line or 'blinded' in line:
        continue
    if 'MET' in line:
      line = "analyze_mutau "+ line + " DATA \n"
    else:
      line = "analyze_mutau "+ line + " MC \n"
  
    outFile.write(line)
   
outFile.close()

# crate signal list in sample_signal.txt
outFile = open("samples_signal.txt", "w")
samples = sorted(os.listdir('/hdfs/store/user/jmadhusu/with_boostedtau/2017_skimmed/with_boostedtaus/signal/'))

for line in samples:
    if '.root' not in line[-6:]:
      continue
    line = line.strip().replace('.root', '')
    line = "analyze_mutau "+ line + " MC \n"
  
    outFile.write(line)
     
outFile.close()
print ("""

samples.txt created
samples_signal.txt created

""")
