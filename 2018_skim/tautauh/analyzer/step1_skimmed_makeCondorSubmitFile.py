import os

inputDrMC="/hdfs/store/user/jmadhusu/2018_skimmed/mutau/"
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
    if "SingleMuon" in j:
        continue
    outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+j+" -1 1000 2018_test MC "+j[:-5]+" $outDir"+"\n")
    if "DYJetsToLL_00" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"DYJetsToLL_Incl_HT_00.root"+" -1 1000 2018_test MC "+"DYJetsToLL_Incl_HT_00"+" $outDir"+"\n")
    if "DYJetsToLL_01" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"DYJetsToLL_Incl_HT_01.root"+" -1 1000 2018_test MC "+"DYJetsToLL_Incl_HT_01"+" $outDir"+"\n")
    if "WJetsToLNu_Incl_08" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"WJetsToLNu_Incl_HT_08.root"+" -1 1000 2018_test MC "+"WJetsToLNu_Incl_HT_08"+" $outDir"+"\n")
    if "WJetsToLNu_Incl_07" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"WJetsToLNu_Incl_HT_07.root"+" -1 1000 2018_test MC "+"WJetsToLNu_Incl_HT_07"+" $outDir"+"\n")
    if "WJetsToLNu_Incl_06" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"WJetsToLNu_Incl_HT_06.root"+" -1 1000 2018_test MC "+"WJetsToLNu_Incl_HT_06"+" $outDir"+"\n")
    if "WJetsToLNu_Incl_05" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"WJetsToLNu_Incl_HT_05.root"+" -1 1000 2018_test MC "+"WJetsToLNu_Incl_HT_05"+" $outDir"+"\n")
    if "WJetsToLNu_Incl_04" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"WJetsToLNu_Incl_HT_04.root"+" -1 1000 2018_test MC "+"WJetsToLNu_Incl_HT_04"+" $outDir"+"\n")
    if "WJetsToLNu_Incl_03" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"WJetsToLNu_Incl_HT_03.root"+" -1 1000 2018_test MC "+"WJetsToLNu_Incl_HT_03"+" $outDir"+"\n")
    if "WJetsToLNu_Incl_02" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"WJetsToLNu_Incl_HT_02.root"+" -1 1000 2018_test MC "+"WJetsToLNu_Incl_HT_02"+" $outDir"+"\n")
    if "WJetsToLNu_Incl_01" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"WJetsToLNu_Incl_HT_01.root"+" -1 1000 2018_test MC "+"WJetsToLNu_Incl_HT_01"+" $outDir"+"\n")
    if "WJetsToLNu_Incl_00" in j :
        outFile.write("./MakeCondorFiles.csh analyze_mutau "+inputDrMC+j+" "+"WJetsToLNu_Incl_HT_00.root"+" -1 1000 2018_test MC "+"WJetsToLNu_Incl_HT_00"+" $outDir"+"\n")


inputDrDATA="/hdfs/store/user/jmadhusu/2018_skimmed/mutau/"
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
check submit.sh
do 

 bash submit.sh


""")

outFile.close()


