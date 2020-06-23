import os
import sys
import fileinput 


with open('out_.txt','w') as write_obj:
    
    for line in fileinput.input('sample_Name.txt'):
        line = line.split(' ')
        rstr = "'"+line[0]+"' : ['"+line[1]+"', '"+line[2]+"', '"+line[3]+"'], \n"
        write_obj.write(rstr)
