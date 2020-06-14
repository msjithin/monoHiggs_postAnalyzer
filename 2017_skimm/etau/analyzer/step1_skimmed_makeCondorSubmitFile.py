import os

inputDrMC="/hdfs/store/user/jmadhusu/2017_skimmed/etau/"
filelistMC=os.listdir(inputDrMC)
filelistMC=sorted(filelistMC)

outFile = open("submit.sh", "w")
outFile.write("""
outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

""")


outFile.write('###########################   MC  #########################'+'\n\n')
outFile.write('./rootcom etau_analyzer analyze_etau  '+'\n\n\n')
for j in filelistMC :
    if 'SingleElectron' in j or '.root' not in j:
        continue
   # outFile.write("./MakeCondorFiles.csh analyze_etau_MC "+inputDrMC+j+" "+"files_initial/"+j+" -1 1000 2018_test MC "+j[:-5]+"\n")
    outFile.write("./MakeCondorFiles.csh analyze_etau "+inputDrMC+j+" "+j+" -1 1000 2017 MC "+j[:-5]+" $outDir"+"\n")
        

inputDrDATA="/hdfs/store/user/jmadhusu/2017_skimmed/etau/"
filelistDATA=os.listdir(inputDrDATA)
filelistDATA=sorted(filelistDATA)

outFile.write('\n\n')
outFile.write('###########################  DATA #########################'+'\n\n')
outFile.write('./rootcom etau_analyzer analyze_etau  '+'\n\n\n')
for j in filelistDATA :
    if 'SingleElectron' in j:
        #outFile.write("./MakeCondorFiles.csh analyze_etau_data "+inputDrDATA+j+" "+"files_initial/"+j+" -1 1000 2018_test DATA "+j[:-5]+"\n")
        outFile.write("./MakeCondorFiles.csh analyze_etau "+inputDrDATA+j+" "+j+" -1 1000 2017 DATA "+j[:-5]+" $outDir"+"\n")
                
print("""
check submit.sh
do 

 bash submit.sh


""")

outFile.close()


