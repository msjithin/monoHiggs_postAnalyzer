from os import listdir

parent_dir = '/hdfs/store/user/jmadhusu/MC2017_12Apr2018_monoHiggs_27Aug2022/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_DYJetsToLL_M-50_TuneCP5_v1/220827_053032'


print listdir(parent_dir)

f = open('condorsubmit.sh', 'w')

f.write("""
./rootcom mutau_analyzer analyze_mutau
outDir="Out_Signal_$(date +"%d-%m-%Y_%H-%M")" 
mkdir $outDir 



""")

fcount = 0
for subdir in listdir(parent_dir):
    for file in listdir(parent_dir+'/'+subdir):
        file_path = parent_dir.replace('/hdfs', 'root://cmsxrootd.hep.wisc.edu/')+'/'+subdir+'/'+file
        index = subdir +'_'+ file.replace('.root', '').split('_')[-1]
        text_to_write = './MakeCondorFiles.csh analyze_mutau {a} DYJetsToLL_M-50_TuneCP5_v1_{b}.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5_v1_{b} $outDir \n'
        f.write( text_to_write.format(a=file_path, b=index) )
        fcount += 1

f.write("""

"""
)

f.close()

print 'total files to be expected = ' , fcount

