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
if histoname=="elept":
  column=1
if histoname=="eleeta":
  column=2
if histoname=="elephi":
  column=3
if histoname=="taupt":
  column=4
if histoname=="taueta":
  column=5
if histoname=="tauphi":
  column=6
if histoname=="eleid":
  column=7
if histoname=="elereliso":
  column=8
if histoname=="tauid":
  column=9
if histoname=="tautightVSe":
  column=10
if histoname=="tauVlooseVSmu":
  column=11
if histoname=="elecharge":
  column=12
if histoname=="taucharge":
  column=13
if histoname=="deltaR":
  column=14
if histoname=="genmatch":
  column=15

print (histoname, column)
def describe_helper(series):
    splits = str(series.describe()).split()
    keys, values = "", ""
    for i in range(0, len(splits), 2):
        keys += "{:8}\n".format(splits[i])
        values += "{:>8}\n".format(splits[i+1])
    return keys, values

myEvents=[]
myPt=[]
fp_mine=  open('../compare_ptdiff_etau.txt', 'r')
for x in fp_mine:
    split_string = x.split()
    #print(split_string[0])
    myEvents.append(float(split_string[0]))
    myPt.append(float(split_string[column]))

thEvents=[]
thPt=[]
fp_their=  open('../compare_ptdiff_smhtt.txt', 'r')
for x in fp_their:
    split_string = x.split()
    #print(split_string[0])
    thEvents.append(float(split_string[0]))
    thPt.append(float(split_string[column]))
    
    
diff=[]
for i in range(len(thEvents)):
    #print( myEvents[i] , thEvents[i] )
    myindex=-1
    if( thEvents[i] in myEvents ):
        myindex = myEvents.index(thEvents[i])
        diff.append(float(myPt[myindex]) - float(thPt[i]))
    
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
plt.ylim(0, 20)
plt.grid(True)
#plt.show()
plt.savefig('plots/pt_diff_'+histoname+'.png')
