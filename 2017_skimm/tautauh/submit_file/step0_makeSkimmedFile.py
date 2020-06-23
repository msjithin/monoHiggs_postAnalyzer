import os
import re

os.popen('ls -d /hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_26April2020/*/*/*/* > mc_dir_list')
inputFile=open("mc_dir_list", "r")
outFile = open("do_skim_mt_mc.sh", "w")
outFile.write("""
./rootcom skimm_mt_2017 analyze_mutau_skim 
outDir="Out_MC_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

""")

def replaceAll(line, searchExp,replaceExp):
    if searchExp in line:
        line = line.replace(searchExp,replaceExp)
    return line
def replaceEnd(line):
    line=line.strip()
    searchExp=line[-6:]
    sample = re.search('/crab_job_(.*)/20', line)
    sampleName = sample.group(1)
    replaceExp = searchExp+"/ "+sampleName+"_"+searchExp[-2:]+".root -1 1000 2017 MC "+sampleName+"_"+searchExp[-2:]+" $outDir"+" \n"
    line = line.replace(searchExp,replaceExp)
    return line

for line in inputFile:
    line =replaceAll(line, "/hdfs/store/","./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/")
    line = replaceEnd(line)
    outFile.write(line)



inputFile.close()
outFile.close()
print """
created do_skim_mt_mc.sh
do 
  bash do_skim_mt_mc.sh
in /nfs_scratch/ directory with required files
"""


os.popen('ls -d /hdfs/store/user/jmadhusu/data2017_31Mar2018_26April2020/SingleMuon/*/*/* > data_dir_list')
inputFile=open("data_dir_list", "r")
outFile = open("do_skim_mt_data.sh", "w")
outFile.write("""
./rootcom skimm_mt_2017 analyze_mutau_skim 
outDir="Out_DATA_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

""")

def replaceAll(line, searchExp,replaceExp):
    if searchExp in line:
        line = line.replace(searchExp,replaceExp)
    return line
def replaceEnd(line):
    line=line.strip()
    searchExp=line[-6:]
    sample = re.search('/crab_job_(.*)/20', line)
    sampleName = sample.group(1)
    replaceExp = searchExp+"/ "+sampleName+"_"+searchExp[-2:]+".root -1 1000 2017 DATA "+sampleName+"_"+searchExp[-2:]+" $outDir"+" \n"
    line = line.replace(searchExp,replaceExp)
    return line

for line in inputFile:
    line =replaceAll(line, "/hdfs/store/","./MakeCondorFiles.csh analyze_mutau_skim /hdfs/store/")
    line = replaceEnd(line)
    outFile.write(line)



inputFile.close()
outFile.close()
print """
created do_skim_mt_data.sh
do 
  bash do_skim_mt_data.sh
in /nfs_scratch/ directory with required files
"""
