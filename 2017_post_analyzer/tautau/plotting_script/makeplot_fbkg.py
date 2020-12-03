#!/usr/bin/env python
import ROOT
import re
from array import array
from sys import argv
import ROOT as Root
import csv
from math import sqrt
from math import pi
import datetime
import argparse
from sys import path
path.append("../../../MacrosAndScripts/")

from myPlotStyle import *


parser = argparse.ArgumentParser()
parser.add_argument("-name",
                    help="name of hist to be plotted")
parser.add_argument("-year",
                    help="dataset year")
parser.add_argument("-in", "--inputFile",
                    help="name of input file")

parser.add_argument("-xaxis",
                    help="x axis label")

parser.add_argument("-lY", "--logYaxis", action="store_true",
                    help="set y axis log scale")

parser.add_argument("-pi", "--piAxis", action="store_true",
                    help="set x axis in intervals of pi")

parser.add_argument("-nL", "--noLegend", action="store_true",
                    help="remove legend from histogram")

parser.add_argument("-cat", "--category",
                    help="pick reco histogram category")

parser.add_argument("-ch", "--channel",
                    help="pick decay channel")

now = datetime.datetime.now()
args =  parser.parse_args()
histoname=args.name
year_=args.year
inputFile_=args.inputFile
category_=args.category
xaxis_label=args.xaxis[:-2]
channel_=args.channel
eventsPerBin=""
if ("muPt" in histoname or "tauPt" in histoname or "elePt" in histoname):
  eventsPerBin=" /2GeV"
if ("higgsPt" in histoname):
  eventsPerBin=" /10GeV"
if ("visMass" in  histoname):
  eventsPerBin=" /5GeV"
if ("met" in histoname):
  eventsPerBin=" /10GeV"
if ("mT_" in histoname):
  eventsPerBin=" /10GeV"
#print "Item / histogram name = ", histoname

if (args.piAxis):
  piAxis = 1
else :
  piAxis = 0
if ( "Phi" in histoname) :
  piAxis = 1

yaxisLog = 0
if (args.logYaxis):
  yaxisLog = 1

noLegend = 0
if (args.noLegend):
  noLegend = 1

 
ROOT.gStyle.SetFrameLineWidth(1)
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetOptStat(0)

c.cd()

if channel_=="mutau":
  OutFile=ROOT.TFile(inputFile_,"r")
if channel_=="etau":
  OutFile=ROOT.TFile(inputFile_,"r") 
if channel_=="tautau":
  OutFile=ROOT.TFile(inputFile_,"r")
if channel_=="emu":
  OutFile=ROOT.TFile(inputFile_,"r") 
if channel_=="combined":
  OutFile=ROOT.TFile(inputFile_,"r") 

adapt=ROOT.gROOT.GetColor(12)
new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)
dirname = [histoname+"/",histoname, histoname+"_fr/", histoname+"_dyll/", histoname+"_dyll_fr/"]

#if (OutFile.GetListOfKeys().Contains(histoname)):
#  print(histoname+' histogram found')
#   OutFile.cd(dirname[1])
#   if (ROOT.gDirectory.GetListOfKeys().Contains("ZTT_"+histoname)):
#     print(dirname[0]+"ZTT_"+histoname+' histogram found')
#   else:
#     print(dirname[0]+"ZTT_"+histoname+' histogram **********NOT********** found')
# else:
#   print(histoname+' histogram **********NOT********** found')
#if (OutFile.Get(dirname[0]+"data_obs_"+histoname)):
#  print(histoname+' histogram found')
#ZTTselect="ZTT"
ZTTselect="ZTTjet"
#WJetselect="WJets"
WJetselect="WJets_jets"
OutFile.cd()
Data_hist    =OutFile.Get(dirname[2]+"data_obs_"+histoname+"_fr")
ZTT_hist     = OutFile.Get(dirname[2]+ZTTselect+"_"+histoname+"_fr")
#ZLL_hist     =  OutFile.Get(dirname[4]+ZTTselect+"_"+histoname+"_dyll_fr")
EWKWMinus_hist = OutFile.Get(dirname[2]+"EWKWMinus_"+histoname+"_fr")
EWKWPlus_hist =  OutFile.Get(dirname[2]+"EWKWPlus_"+histoname+"_fr")
EWKZ2Jets_hist = OutFile.Get(dirname[2]+"EWKZ2Jets_"+histoname+"_fr")
#Wjets_hist   = 
GluGluH_hist = OutFile.Get(dirname[2]+"GluGluH_"+histoname+"_fr")
ST_t_hist    = OutFile.Get(dirname[2]+"ST_t_"+histoname+"_fr")
TT_hist      = OutFile.Get(dirname[2]+"TT_"+histoname+"_fr")
VBFH_hist    =  OutFile.Get(dirname[2]+"VBFH_"+histoname+"_fr")
VV_hist      =  OutFile.Get(dirname[2]+"VV_"+histoname+"_fr")
WminusH_hist = OutFile.Get(dirname[2]+"WminusH_"+histoname+"_fr")
WplusH_hist  = OutFile.Get(dirname[2]+"WplusH_"+histoname+"_fr")
ZH_hist      = OutFile.Get(dirname[2]+"ZH_"+histoname+"_fr")
VVV_hist     = OutFile.Get(dirname[2]+"VVV_"+histoname+"_fr")
#ZJetsToNuNu_hist = 


if( OutFile.Get(dirname[0]+"EWKWMinus_"+histoname) ) :
  GluGluH_hist.Add(EWKWMinus_hist)
if( OutFile.Get(dirname[0]+"EWKWPlus_"+histoname) ):
  GluGluH_hist.Add(EWKWPlus_hist)
if( OutFile.Get(dirname[0]+"EWKZ2Jets_"+histoname) ):
  GluGluH_hist.Add(EWKZ2Jets_hist)
if(  OutFile.Get(dirname[0]+"VBFH_"+histoname) ):
  GluGluH_hist.Add(VBFH_hist)
if(OutFile.Get(dirname[0]+"WminusH_"+histoname)):
  GluGluH_hist.Add(WminusH_hist)
if(OutFile.Get(dirname[0]+"WplusH_"+histoname)):  
  GluGluH_hist.Add(WplusH_hist)
if( OutFile.Get(dirname[0]+"ZH_"+histoname)):
  GluGluH_hist.Add(ZH_hist)
if(OutFile.Get(dirname[0]+"ST_t_"+histoname)):
  GluGluH_hist.Add(ST_t_hist)
if( OutFile.Get(dirname[0]+"VVV_"+histoname) ):
  GluGluH_hist.Add(VVV_hist)
if(OutFile.Get(dirname[0]+"VV_"+histoname)):
  GluGluH_hist.Add(VV_hist)


#sampleList    = [Data_hist,    ZTT_hist,   TT_hist,   GluGluH_hist,     ]
#sampleListRef = ['Data_hist', 'ZTT_hist',  'TT_hist', 'GluGluH_hist'  ]
sampleList    = [Data_hist,   GluGluH_hist  , ZTT_hist  ]
sampleListRef = ['Data_hist', 'GluGluH_hist','ZTT_hist' ]
if channel_=="etau":
  sampleList    = [Data_hist,    ZTT_hist,   TT_hist,  GluGluH_hist ]
  sampleListRef = ['Data_hist', 'ZTT_hist',  'TT_hist', 'GluGluH_hist']
# if channel_=="etau":x1
#  sampleList    = [Data_hist,    ZTT_hist,  ZLL_hist, TT_hist,   GluGluH_hist ]
#  sampleListRef = ['Data_hist', 'ZTT_hist', 'ZLL_hist', 'TT_hist', 'GluGluH_hist']


#Wjets_hist.SetFillColor(ROOT.TColor.GetColor(color_wjets))
#F_bkg.SetFillColor(ROOT.TColor.GetColor(color_wjets))
ZTT_hist.SetFillColor(ROOT.TColor.GetColor(color_ztt))
#ZLL_hist.SetFillColor(ROOT.TColor.GetColor(color_zll))
TT_hist.SetFillColor(ROOT.TColor.GetColor(color_tt))
GluGluH_hist.SetFillColor(ROOT.TColor.GetColor(color_ggh))
VV_hist.SetFillColor(ROOT.TColor.GetColor(color_vv))

for i in range(len(sampleList)):
  sampleList[i].SetLineColor(1)

stack=ROOT.THStack("stack","stack")  
for i in range(len(sampleList)-1, 0, -1):
  #print(i, sampleList[i])
  stack.Add(sampleList[i])

errorBand=sampleList[1].Clone()
for i in range(2, len(sampleList)):
  errorBand.Add(sampleList[i])

#errorBand.SetMarkerSize(0)
errorBand.SetFillColor(1)
errorBand.SetFillStyle(errorStyle)
#errorBand.SetLineWidth(1)

Data_hist.GetXaxis().SetTitle("")
Data_hist.GetXaxis().SetTitleSize(0)

nDivXAxis=405
if piAxis == 1:
  nDivXAxis= -105
elif histoname=='cutflow_n' or histoname=='cutflow_Htt':
  nDivXAxis= Data_hist.GetNbinsX()
elif histoname=='relMuIso_5' :
  nDivXAxis= Data_hist.GetNbinsX()
elif histoname=='tauAntiEle_5' or histoname=='tauAntiMu_5':
  nDivXAxis= Data_hist.GetNbinsX()
elif histoname=='muCharge_5' or histoname=='tauCharge_5':
  nDivXAxis= Data_hist.GetNbinsX()
elif histoname=='tauDecayMode_5':
  nDivXAxis= Data_hist.GetNbinsX()

Data_hist.GetXaxis().SetNdivisions(nDivXAxis)
Data_hist.GetXaxis().SetLabelSize(0)
Data_hist.GetYaxis().SetLabelFont(42)
Data_hist.GetYaxis().SetLabelOffset(0.01)
Data_hist.GetYaxis().SetLabelSize(0.04)
Data_hist.GetYaxis().SetTitleSize(0.05)
Data_hist.GetYaxis().SetTitleOffset(1.22)
Data_hist.SetTitle("")
Data_hist.GetYaxis().SetTitle("")


c.cd()
pad1.Draw()
pad1.cd()

if yaxisLog == 1 :
  pad1.SetLogy()  

Data_hist.SetMarkerStyle(20)
Data_hist.SetMarkerColor(1)
Data_hist.SetMarkerSize(1.5)
if histoname == "cutflow_n":
  Data_hist.SetMarkerSize(2)


Data_hist.GetXaxis().SetTitle("")
Data_hist.GetYaxis().SetTitle("Events"+eventsPerBin)
if yaxisLog == 1 :
  Data_hist.SetMaximum(100*max(Data_hist.GetMaximum(),stack.GetMaximum()))
  Data_hist.SetMinimum(1000)
else :
  if channel_=="mutau":
    Data_hist.SetMaximum(1.055*max(Data_hist.GetMaximum(),stack.GetMaximum()))
  if channel_=="etau":
    Data_hist.SetMaximum(1.5*max(Data_hist.GetMaximum(), stack.GetMaximum()))
  if channel_=="tautau":
    Data_hist.SetMaximum(1.55*max(Data_hist.GetMaximum(), stack.GetMaximum()))
  Data_hist.SetMinimum(0.0)
  
if histoname == "cutflow_n":
  Data_hist.SetMinimum(1000)

Data_hist.Draw("e1")
stack.Draw("histsame")
# if histoname != "cutflow_n":
#   errorBand.Draw("e1same")
Data_hist.Draw("e1same")
if histoname == "cutflow_n":
  Data_hist.Draw("e0psame")

legendNameList = {
  'Data_hist'  : 'Data obs',
  'ZTT_hist'   : 'Z->tautau',
  'ZLL_hist'   : 'Z-> ll',
  'Wjets_hist' : 'WJets',
  'F_bkg'      : 'jet-tau fake', 
  'TT_hist'    : 'ttabr',
  #'GluGluH_hist' : 'ggh, vbfH, ZH',
  'GluGluH_hist' : 'Others',
  'VV_hist'    : 'VV, SingleTop',
  'ZJetsToNuNu_hist' : 'Z->nunu + jets'
}
legende=make_legend()
for i in range(len(sampleListRef)):
  if(i==0):
    legende.AddEntry(sampleList[i], legendNameList[sampleListRef[i]], "elp")
  else:
    legende.AddEntry(sampleList[i], legendNameList[sampleListRef[i]], "f")
  

l1=add_lumi(year_, channel_)
l1.Draw("same")
l2=add_CMS()
l2.Draw("same")
l3=add_Preliminary()
l3.Draw("same")

pad1.RedrawAxis()
#print "Line 217 is okay"


c.cd()
pad2.Draw()
pad2.cd()
h1=Data_hist.Clone()
h3=errorBand.Clone()
hwoE=errorBand.Clone()
h1.SetMarkerStyle(20)
h1.SetMarkerSize(2.0)
#h1.Sumw2()
#h3.Sumw2()
h1.SetStats(0)
h1.Divide(hwoE)
h3.Divide(hwoE)
h1.GetXaxis().SetTitle(xaxis_label)
h1.SetMarkerColor(1)
h1.SetLineColor(1)
h1.SetTitle("")
h1.GetXaxis().SetLabelSize(0.1)
h1.GetYaxis().SetTitle("Data/MC")
h1.GetYaxis().SetLabelSize(0.08)
h1.SetMaximum(1.5)
h1.SetMinimum(0.5)
h1.GetXaxis().SetNdivisions(nDivXAxis)
if piAxis ==1 :
  h1.GetXaxis().SetBinLabel(30,"#pi");
  h1.GetXaxis().SetBinLabel(22,"#frac{#pi}{2}");
  h1.GetXaxis().SetBinLabel(15,"0");
  h1.GetXaxis().SetBinLabel(8,"#frac{-#pi}{2}");
  h1.GetXaxis().SetBinLabel(1,"-#pi");
  
if piAxis == 1:
  h1.GetXaxis().SetLabelOffset(0.1)
  h1.GetXaxis().SetLabelSize(0.15)
  h1.GetXaxis().SetNdivisions(-105)

elif histoname=='Events_level_':
  
  h1.GetXaxis().SetBinLabel(13,"METcut");
  h1.GetXaxis().SetBinLabel(12,"VisibleMass");
  h1.GetXaxis().SetBinLabel(11,"HiggsPt");
  h1.GetXaxis().SetBinLabel(10,"Mt cut");
  h1.GetXaxis().SetBinLabel(9,"BjetVeto");
  h1.GetXaxis().SetBinLabel(8,"ThirdLepVeto");
  h1.GetXaxis().SetBinLabel(7,"DeltaR");
  h1.GetXaxis().SetBinLabel(6,"opp. charge");
  h1.GetXaxis().SetBinLabel(5,"GoodTau");
  h1.GetXaxis().SetBinLabel(4,"GoodMuon");
  h1.GetXaxis().SetBinLabel(3,"mu Trg");
  h1.GetXaxis().SetBinLabel(2,"MET Filters");
  h1.GetXaxis().SetBinLabel(1,"No cuts");
  h1.GetXaxis().SetNdivisions(-115)
  h1.GetXaxis().LabelsOption("v")
  h1.GetXaxis().SetTitle(" cut flow")
  h1.SetMaximum(1.2)#FIXME(1.2)
  h1.SetMinimum(0.8)#FIXME(0.8)

elif histoname=='cutflow_n':
  
  h1.GetXaxis().SetBinLabel(8,"mu-tau \n separation");
  h1.GetXaxis().SetBinLabel(7,"b-Jet veto");
  h1.GetXaxis().SetBinLabel(6,"3rd lepton \n veto");
  h1.GetXaxis().SetBinLabel(5,"opposite \n ch.");
  h1.GetXaxis().SetBinLabel(4,"tau \n selections");
  h1.GetXaxis().SetBinLabel(3,"Muon \n selections");
  h1.GetXaxis().SetBinLabel(2,"Trigger");
  h1.GetXaxis().SetBinLabel(1,"Initial");
  h1.GetXaxis().LabelsOption("v")
  h1.GetXaxis().SetTitle(" ")
  h1.SetMaximum(2.0)#FIXME(1.2)
  h1.SetMinimum(0.0)#FIXME(0.8)



h1.GetYaxis().SetNdivisions(5)
h1.GetXaxis().SetTitleFont(42)
h1.GetYaxis().SetTitleFont(42)
h1.GetXaxis().SetTitleSize(0.15)
h1.GetYaxis().SetTitleSize(0.15)
h1.GetYaxis().SetTitleOffset(0.3)
h1.GetXaxis().SetTitleOffset(1.1)

h1.Draw("e0p")
h3.Draw("e2same")

c.cd()
#pad1.RedrawAxis()
#pad2.RedrawAxis()
if noLegend == 0:
  legende.Draw()
#ROOT.gPad.RedrawAxis()

c.Modified()
c.SaveAs("plots/plot_fbkg_"+histoname+"_"+channel_+".png")

