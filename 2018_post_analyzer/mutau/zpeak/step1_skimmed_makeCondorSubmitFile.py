import os

inputDrMC="/hdfs/store/user/jmadhusu/2018_skimmed/muSelections/"
filelistMC=os.listdir(inputDrMC)
filelistMC=sorted(filelistMC)

outFile = open("testAn.sh", "w")
outFile.write("""
outDir="Out_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 

""")


outFile.write('###########################   MC  #########################'+'\n\n')
outFile.write('./rootcom mutau_analyzer analyze_mutau  '+'\n\n\n')
for j in filelistMC :
    if not 'SingleMuon' in j:
        # outFile.write("./MakeCondorFiles.csh analyze_mutau_MC "+inputDrMC+j+" "+"files_initial/"+j+" -1 1000 2018_test MC "+j[:-5]+"\n")
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+j+" -1 1000 2018_test MC "+j[:-5]+" $outDir"+"\n")
        

inputDrDATA="/hdfs/store/user/jmadhusu/2018_skimmed/muSelections/"
filelistDATA=os.listdir(inputDrDATA)
filelistDATA=sorted(filelistDATA)

outFile.write('\n\n')
outFile.write('###########################  DATA #########################'+'\n\n')
outFile.write('./rootcom mutau_analyzer analyze_mutau  '+'\n\n\n')
for j in filelistDATA :
    if 'SingleMuon' in j:
        #outFile.write("./MakeCondorFiles.csh analyze_mutau_data "+inputDrDATA+j+" "+"files_initial/"+j+" -1 1000 2018_test DATA "+j[:-5]+"\n")
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrDATA+j+" "+j+" -1 1000 2018_test DATA "+j[:-5]+" $outDir"+"\n")
                
print("""
check testAn.sh
do 

 bash testAn.sh


""")

outFile.close()


