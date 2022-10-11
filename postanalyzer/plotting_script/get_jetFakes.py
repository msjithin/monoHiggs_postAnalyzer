#!/usr/bin/env python
import ROOT
import re
from array import array
from os import path
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
#mc_samples = ['ZTTjet', 'EWKWMinus', 'EWKWPlus', 'EWKZ2Jets', 'GluGluH', 'GluGluZH', 'HWminusJ', 'HWplusJ', 'HZJ', 'ST_t', 'TT', 'VBFH', 'WGToLNuG', 'VV', 'VVV', 'WplusH', 'ZH', 'ZJetsToNuNu']
mc_samples = ['ZTTjet', 'ZLLjet', 'TT', 'otherMC', 'STT', 'VVT']

        
def checkHistogram(f, histogram):
    isthere=  f.GetListOfKeys().Contains(histogram)
    #print(isthere)
    return isthere


def GetKeyNames( self, dir = "" ):
    self.cd(dir)
    return [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
ROOT.TFile.GetKeyNames = GetKeyNames

def getHistList(inFile, isblinded=False):
    keyList = inFile.GetKeyNames()
    #print "\nKeys in file:", keyList
    tmpList= []
    for tdir in keyList:
        if "_fr" not in tdir: continue
        if "CMS_htt_boson" in tdir : continue
        print("\n "+tdir)

        data_dir  = tdir+'/'+'data_obs_'+tdir
        if isblinded:
            data_dir  = tdir+'/blinded_data_obs_'+tdir
        print(" checking "+data_dir)
        if not inFile.Get(data_dir):
            print("################## No histogram for {} ".format(data_dir))
            print("################## check this please")
            continue
        jetFakes = inFile.Get(data_dir)
        print(tdir+' integral  data_obs', jetFakes.Integral())
        for mc in mc_samples:
            tmppath = tdir+'/'+mc+'_'+tdir
            try:
                tmpHist = inFile.Get(tmppath)
                jetFakes.Add(tmpHist, -1)
            except:
                pass
        print('integral  jetFakes', jetFakes.Integral())
        tdir = tdir.replace('_fr', '')
        inFile.cd(tdir)
        jetFakes.SetName("jetFakes_"+tdir)
        jetFakes.Write()
        
    return 


def main(histogram, isblinded):
    if not path.exists("f_"+histogram+"_initial.root"):
        print("File doesnt exist for {} ".format(histogram))
        exit()

    inFile_nominal= ROOT.TFile("f_"+histogram+"_initial.root","UPDATE")
    getHistList(inFile_nominal, isblinded)
    inFile_nominal.Close()
    

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist",
                    help="name of histogram elePt , tauPt, ..  Default=etau")
    parser.add_argument("--blinded",
                    help="is this for blinded case,   Default= 0",
                    choices=('0', '1'),
                    default='0'
                    )
    args =  parser.parse_args()
    if args.hist is None:
        print("Specify histogram using --hist")
        exit()   
    else:
        histogram = args.hist
    if args.blinded:
        print("running blinded with 1/5th data")
    else:
        print("running with full data")
    print('histogram = ' , histogram)
    isBlinded = False
    if args.blinded=='1' :
        isBlinded = True
    main(histogram, isBlinded)


