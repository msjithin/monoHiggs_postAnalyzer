import os

inputDrMC="/hdfs/store/user/jmadhusu/2017_skimmed/htt_et_2017/"
filelistMC=os.listdir(inputDrMC)
filelistMC=sorted(filelistMC)

outFile = open("submit_etau.sh", "w")
outFile.write("""
outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

""")


outFile.write('###########################   MC  #########################'+'\n\n')
outFile.write('./rootcom etau_analyzer executable_etau  '+'\n\n\n')
for j in filelistMC :
    if 'SingleElectron' in j or '.root' not in j:
        continue
    outFile.write("./MakeCondorFiles.csh executable_etau "+inputDrMC+j+" "+j+" -1 1000 2017 MC "+j[:-5]+" $outDir"+"\n")
    #outFile.write("./executable_etau "+inputDrMC+j+" "+j+" -1 1000 2017 MC "+j[:-5]+" $outDir"+"\n")
    

outFile.close()


