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
#mc_samples = ['ZTTjet', 'EWKWMinus', 'EWKWPlus', 'EWKZ2Jets', 'GluGluH', 'GluGluZH', 'HWminusJ', 'HWplusJ', 'HZJ', 'ST_t', 'TT', 'VBFH', 'WGToLNuG', 'VV', 'VVV', 'WplusH', 'ZH', 'ZJetsToNuNu']


def drawHist(hist, category):
    c=ROOT.TCanvas("canvas","",0,0,1300,1200)
    c.cd()
    try:
        hist.Draw()
        c.SaveAs("vismass_"+category+".png")
    except:
        print("troubel making "+"vismass_"+category+".png" )
    c.Clear()
    c.Close()

        
# #tot_TMass_9_dyll_JER_down
# h_selected = "tot_TMass_full"
# #HistSelected = OutFile.Get("tot_TMass")
def checkHistogram(f, histogram):
    isthere=  f.GetListOfKeys().Contains(histogram)
    #print(isthere)
    return isthere


def GetKeyNames( self, dir = "" ):
    self.cd(dir)
    return [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
ROOT.TFile.GetKeyNames = GetKeyNames

def getHistList(inFile):
    keyList = inFile.GetKeyNames()
    #print "\nKeys in file:", keyList
    tmpList= []
    for tdir in keyList:
        if "_dyll" not in tdir: continue
        
        try:
            ZLL = inFile.Get(tdir+'/'+'ZTTjet_'+tdir)
            print('integral  ZLL', ZLL.Integral())
            #path = path.replace("_dyll", "")
            tdir = tdir.replace("_dyll", "")
            inFile.cd(tdir)
            ZLL.SetName("ZLLjet_"+tdir)
            ZLL.Write()
        except:
            pass
    return 




def main(histogram):
    inFile_nominal= ROOT.TFile("f_"+histogram+"_initial.root","UPDATE")
    inFile_up     = ROOT.TFile("f_"+histogram+"_up.root","UPDATE")
    inFile_down   = ROOT.TFile("f_"+histogram+"_down.root","UPDATE")
    getHistList(inFile_nominal)
    inFile_nominal.Close()
    getHistList(inFile_up)
    inFile_up.Close()
    getHistList(inFile_down)
    inFile_down.Close()


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist",
                    help="name of histogram elePt , tauPt, ..  Default=etau")
    args =  parser.parse_args()
    if args.hist is None:
        histogram = 'etau'
    else:
        histogram = args.hist
    print('histogram = ' , histogram)
    main(histogram)

