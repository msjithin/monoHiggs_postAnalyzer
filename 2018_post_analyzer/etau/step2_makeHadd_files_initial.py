import os

filelist=os.listdir("rootFiles")
samplelist=[]

res = [] 
for x in filelist:
    x=x[:-8]
    if x not in res:
        res.append(x)
    
outFile = open("hadd_files.sh", "w")
outFile.write(
"""
#!/bin/bash
set -e 

if [ "$(ls -A files_initial)" ]; then
echo "Take action files_initial/ is not Empty .... removing existing files ....."
rm files_initial/*.root
else
echo " files_initial/ is Empty"
fi

"""
)

for i in res :
    filestr=''
    for j in filelist :
        
        tmp = j[:-8]
        if len(i)==len(tmp) and (i in j) :
            filestr+='rootFiles/'+j+' '
    outFile.write('hadd files_initial/{}_final.root {}\n'.format(i, filestr))

outFile.write(
"""

"""
)
print("""
output written to hadd_files.sh 
do  ' bash hadd_files.sh '
"""
)
outFile.close()
