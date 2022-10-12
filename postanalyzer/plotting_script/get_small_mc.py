#!/usr/bin/env python
import ROOT
import re
from array import array
import sys
import csv
from math import sqrt
from math import pi
import datetime
import argparse
# from sys import path
# path.append("../../../MacrosAndScripts/")
# from myPlotStyle import *
ROOT.gStyle.SetFrameLineWidth(1)
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetOptStat(0)
#mc_samples = [ 'EWKWMinus', 'EWKWPlus', 'EWKZ2Jets', 'GluGluH', 'GluGluZH', 'ST_t', 'VBFH', 'VV', 'VVV', 'WplusH', 'WminusH', 'ZH', 'ZJetsToNuNu']
mc_dict = {'STT': ['ST_t'] ,
           'VVT': ['VV', 'VVV' ],
           'otherMC' : [ 'EWKWMinus', 'EWKWPlus', 'EWKZ2Jets', 'GluGluH', 'GluGluZH' , 'VBFH', 'WplusH', 'WminusH', 'ZH', 'ZJetsToNuNu' ]
       }
final_mc_list=['ZTTjet', 'ZLLjet', 'TT']
        
def checkHistogram(f, histogram):
    isthere=  f.GetListOfKeys().Contains(histogram)
    #print(isthere)
    return isthere


def GetKeyNames( self, dir = "" ):
    self.cd(dir)
    return [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
ROOT.TFile.GetKeyNames = GetKeyNames


def getHistList_v3(inFile):
    keyList = inFile.GetKeyNames()
    tmpList= []
    for tdir in keyList:
        print('\n\ntdir  =  ', tdir)
        for key, value in mc_dict.items():
            #print key, value
            mc_samples = value
            i = 0
            while i<len(mc_samples):
                if inFile.Get(tdir+'/'+mc_samples[i]+'_'+tdir):
                    break
                i += 1
            if i==len(mc_samples):
                print("tdir = "+tdir+" didnt have any matching cases ")
                print("##################### nothing filled for ", mc_samples)
                continue
            smallMC = inFile.Get(tdir+'/'+mc_samples[i]+'_'+tdir).Clone()
            smallMC.Reset("ICES")
            for mc in mc_samples:
                tmppath = tdir+'/'+mc+'_'+tdir
                if inFile.Get(tmppath):
                    tmpHist = inFile.Get(tmppath)
                    smallMC.Add(tmpHist)                
            inFile.cd(tdir)
            smallMC.SetName(key+"_"+tdir)
            print("Small MC "+key+"_"+tdir+" integral = "+str(smallMC.Integral()))
            smallMC.Write()
    return 

    
def main(histogram):
    inFile_nominal= ROOT.TFile("f_"+histogram+"_initial.root","UPDATE")
    inFile_up     = ROOT.TFile("f_"+histogram+"_up.root","UPDATE")
    inFile_down   = ROOT.TFile("f_"+histogram+"_down.root","UPDATE")
    getHistList_v3(inFile_nominal)
    inFile_nominal.Close()
    getHistList_v3(inFile_up)
    inFile_up.Close()
    getHistList_v3(inFile_down)
    inFile_down.Close()


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist",
                    help="name of histogram elePt , tauPt, ..  Default=etau")
    args =  parser.parse_args()
    if args.hist is None:
        print("No aruments passed, which varable? ")
        exit()
    histogram = args.hist
    print('histogram = ' , histogram)
    main(histogram)

