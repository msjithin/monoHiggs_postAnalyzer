import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sys import argv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-name",
                    help="name of hist to be plotted")
args =  parser.parse_args()
histoname=args.name

column=0
if histoname=="event":
    column=0
if histoname=="lumi":
    column=1
if histoname=="run":
    column=2
if histoname=="elept":
    column=3
if histoname=="taupt":
    column=4
if histoname=="met":
    column=5
if histoname=="taudm":
    column=6
if histoname=="tauGM":
    column=7

myFile = open('../eventAnalysis_etau.txt', 'r')
myLines = myFile.readlines()

# thFile =  open('eventAnalysis.txt', 'r')
# thLines = thFile.readlines()
def describe_helper(series):
    splits = str(series.describe()).split()
    keys, values = "", ""
    for i in range(0, len(splits), 2):
        keys += "{:8}\n".format(splits[i])
        values += "{:>8}\n".format(splits[i+1])
    return keys, values

myEvents=[]
thEvents=[]

diff=[]
for line in myLines:
    myList=line.strip().split("\t")
    #print myList
    if  myList[0]=='event':
        continue
    with open('../eventAnalysis.txt') as f:
        for thline in f:
            thList=thline.strip().split("\t")
            #print thList
            if( myList[0]==thList[0]):
                # print "common event :" + myList[0] +" "+ thList[0]
                #print myList
                # print thList
                myEvents.append(float(myList[column]))
                thEvents.append(float(thList[column]))
                diff.append(float(myList[column])-float(thList[column]))
                

# print myEvents
# print thEvents
# print diff_array
diff_array=np.array(diff)

plt.hist(diff, bins=50, )
#plt.figtext(1.0, 0.2, diff_array.describe())
plt.figtext(0.70, .49, describe_helper(pd.Series(diff_array))[0], {'multialignment':'left'})
plt.figtext(0.80, .49, describe_helper(pd.Series(diff_array))[1], {'multialignment':'left'})

title_=histoname
plt.title(title_+' difference ')
plt.ylabel("number of events")
plt.xlabel(histoname+" 1 - "+histoname+ "2")
plt.xlim(-5, 5)
plt.ylim(0, 220)
plt.grid(True)
#plt.show()
plt.savefig('pt_diff_'+histoname+'.png')

myFile.close()
#thFile.close()
