import os
import sys
import fileinput 


with open('out_.txt','w') as write_obj:
    
    for line in fileinput.input('list'):
        line=line.strip('\n')
        sampleName=line.split('_')
        argEnd="0"
        rstr = "'"+line+"' : ['"+sampleName[0]+"', '"+sampleName[0]+"', '"+argEnd+"'], "
        print(rstr)
        #write_obj.write(rstr)
