import ROOT
from array import array
from sys import exit
from os import listdir
import re
import argparse
import lumi_weights_2017 as lumi
import argparse
from main import var_mapping


def GetKeyNames( self, dir = "" ):
    self.cd(dir)
    return [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
ROOT.TFile.GetKeyNames = GetKeyNames

def getSaveName(histName):
    position = re.search(r"_\d", histName)
    saveName = ""
    #print 'histName = ', histName , 'position=', position.start()
    #print histName[: position.start()]
    if position:
        saveName =  histName[: position.start()] if histName[position.start()]=="_" else histName[: position.start()+1]
    else:
        saveName =  histName
    if 'up' in histName:
        saveName += "_up"
    elif 'down' in histName:
        saveName += "_down"
    #print 'saveName = ', saveName
    #print ".................. \n"
    return saveName

def getHistList(sampleName = "", hist_name="", idx=""):
    inFile= ROOT.TFile("../files_initial/"+sampleName,"r")
    keyList = inFile.GetKeyNames()
    print "{} histograms in file {}".format(sampleName, len(keyList))
    hist_mapping = { }

    search_str = hist_name+'_'+idx
    print 'Searching for ',  search_str
    
    for hist in keyList:

        result1 = hist.find(search_str)    
        result2 = hist.find(search_str+'_')
        if not result1+len(search_str)==len(hist) and not result2>=0:       
            continue
        saveName = getSaveName(hist)
        if saveName in hist_mapping:
            hist_mapping[saveName].append(hist)
        else:
            hist_mapping[saveName] = [hist]
    
    inFile.Close()
    if not hist_mapping:
        print('\nrequested histogram name not in file....................')
    # for k, v in hist_mapping.items():
    #     print k, v
    return hist_mapping
    

def make_files(sampleName = "", hist_name="", idx="", isBlinded=False):
    inFile= ROOT.TFile("../files_initial/"+sampleName,"r")
    nEventsHisto = inFile.Get("nEvents")
    if not isinstance(nEventsHisto, ROOT.TH1F):
        print 'nEvents not found in ' "../files_initial/"+sampleName
        return 
    nGeneratedEvents = nEventsHisto.GetBinContent(1)
    weight, saveName= lumi.get_lumiweight(sampleName[:-11], nGeneratedEvents, isBlinded)
    
    hist_mapping =  getHistList(sampleName, hist_name, idx)
    #print ""
    #print hist_mapping.keys()
    #print ""
     
    for key in hist_mapping:
        outFile = ROOT.TFile("sample/"+sampleName.replace('.root', '')+'_'+key+'.root' ,  "UPDATE")
        #print 'in file ', "sample/"+sampleName[:-5]+'_'+key+'.root'
        for histName in hist_mapping[key]:
            outFile.cd()
            tmpHist =  inFile.Get(histName)
            if not outFile.GetDirectory(histName):
                outFile.mkdir(histName)
            outFile.cd(histName)
            tmpHist.Scale(weight)
            tmpHist.Write(saveName+"_"+histName)
        outFile.Close()
    inFile.Close()


def postAnalyzer(hist, idx, isBlinded):
    path = "../files_initial/"
    f_list = sorted(listdir(path))
    
    for infile in f_list:
        #check_integral(infile)
        if '.root' not in infile: continue
        print "Making root files for "+var_mapping[int(hist)]
        make_files(infile , var_mapping[int(hist)], idx, isBlinded)
        print 'Done......................'
    

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist",
                    help="index of hist to be plotted,  example -idx 4   or -idx 4,5,12,11 ")
    parser.add_argument("--idx",
                    help="index of selection to be plotted,  example -idx 9, 9b")
    parser.add_argument("--blinded",
                    help="is this for blinded case,   Default= 0",
                    choices=('0', '1'),
                    default='0'
                    )
    args =  parser.parse_args()
    if args.idx is None or args.hist is None:
        print "No index passed"
        print "USAGE : python postAnalyzer.py -hist 1 -idx 9 "
        exit('No index passed..............')
    hist = args.hist
    idx = args.idx
    isBlinded = False
    if args.blinded=='1' :
        isBlinded = True
    postAnalyzer(hist, idx, isBlinded)
