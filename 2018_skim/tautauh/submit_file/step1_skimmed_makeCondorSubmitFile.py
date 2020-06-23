import os

inputDrMC="/hdfs/store/user/jmadhusu/2017_skimmed/mutau/"
filelistMC=os.listdir(inputDrMC)
filelistMC=sorted(filelistMC)

outFile = open("submit.sh", "w")
outFile.write("""
outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

""")


outFile.write('###########################   MC  #########################'+'\n\n')
outFile.write('./rootcom mutau_analyzer analyze_mutau  '+'\n\n\n')
for j in filelistMC :
    if 'SingleMuon' in j or '.root' not in j:
        continue
   # outFile.write("./MakeCondorFiles.csh analyze_mutau_MC "+inputDrMC+j+" "+"files_initial/"+j+" -1 1000 2018_test MC "+j[:-5]+"\n")
    outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+j+" -1 1000 2017 MC "+j[:-5]+" $outDir"+"\n")
        

inputDrDATA="/hdfs/store/user/jmadhusu/2017_skimmed/mutau/"
filelistDATA=os.listdir(inputDrDATA)
filelistDATA=sorted(filelistDATA)

outFile.write('\n\n')
outFile.write('###########################  DATA #########################'+'\n\n')
outFile.write('./rootcom mutau_analyzer analyze_mutau  '+'\n\n\n')
for j in filelistDATA :
    if 'SingleMuon' in j:
        #outFile.write("./MakeCondorFiles.csh analyze_mutau_data "+inputDrDATA+j+" "+"files_initial/"+j+" -1 1000 2018_test DATA "+j[:-5]+"\n")
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrDATA+j+" "+j+" -1 1000 2017 DATA "+j[:-5]+" $outDir"+"\n")
                
print("""
check submit.sh
do 

 bash submit.sh


""")

outFile.close()


