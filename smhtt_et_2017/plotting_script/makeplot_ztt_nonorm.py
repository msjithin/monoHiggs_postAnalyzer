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
# from sys import path
# path.append("../../../MacrosAndScripts/")

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
parser.add_argument("-bkg", "--background",
                    help="pick MC bkg to be plotted")


now = datetime.datetime.now()
args =  parser.parse_args()
histoname=args.name
year_=args.year
inputFile_=args.inputFile
category_=args.category
xaxis_label=args.xaxis[:-2]
channel_=args.channel
bkg_=args.background
bkg_selected=bkg_
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

# if channel_=="mutau":
#   OutFile=ROOT.TFile(inputFile_,"r")
# if channel_=="etau":
#   OutFile=ROOT.TFile(inputFile_,"r") 
# if channel_=="tautau":
#   OutFile=ROOT.TFile(inputFile_,"r")
# if channel_=="emu":
#   OutFile=ROOT.TFile(inputFile_,"r") 
# if channel_=="combined":
#   OutFile=ROOT.TFile(inputFile_,"r") 
histSelect=bkg_selected
OutFile=ROOT.TFile("../mine_rootfile/"+histSelect+".root","r")
OutFile_1=ROOT.TFile("../theirs_rootfile/"+histSelect+".root","r")




adapt=ROOT.gROOT.GetColor(12)
new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)
dirname = [histoname+"/",histoname, histoname+"_fr/", histoname+"_dyll/", histoname+"_dyll_fr/"]


HistSelected = OutFile.Get(histoname)
HistSelected_1 = OutFile_1.Get(histoname)
Data_hist    = OutFile.Get(histoname)

dy_norm=2.447677534
# if(bkg_=="DY"):
#   HistSelected.Scale(dy_norm)
#   HistSelected_1.Scale(dy_norm)
#   Data_hist.Scale(dy_norm)

HistSelected.SetFillColor(ROOT.TColor.GetColor(color_ztt))
#HistSelected_1.SetFillColor(1)


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

with open('eventYield.csv', mode='a') as yield_file:
  yield_write = csv.writer(yield_file, delimiter=',', quotechar='"')
  yield_write.writerow([ histoname ])
  yield_write.writerow(['monoH', Data_hist.Integral() ])
  yield_write.writerow(['SMHtt_et', HistSelected_1.Integral() ])  
c.cd()
pad1.Draw()
pad1.cd()

if yaxisLog == 1 :
  pad1.SetLogy()  

Data_hist.SetMarkerStyle(20)
Data_hist.SetMarkerSize(0)
HistSelected_1.SetMarkerStyle(20)
HistSelected_1.SetMarkerSize(2)
HistSelected_1.SetMarkerColor(1)

if histoname == "cutflow_n":
  Data_hist.SetMarkerSize(0)
Data_hist.SetMarkerColor(1)

Data_hist.GetXaxis().SetTitle("")
Data_hist.GetYaxis().SetTitle("Events"+eventsPerBin)
if yaxisLog == 1 :
  Data_hist.SetMaximum(100*Data_hist.GetMaximum())
  Data_hist.SetMinimum(1000)
else :
  if channel_=="mutau":
    Data_hist.SetMaximum(1.055*max(Data_hist.GetMaximum(), HistSelected_1.GetMaximum()))
  if channel_=="etau":
    Data_hist.SetMaximum(1.5*max(Data_hist.GetMaximum(), HistSelected_1.GetMaximum()))
  Data_hist.SetMinimum(0.0)
if histoname == "cutflow_n":
  Data_hist.SetMinimum(1000)

errorBand=HistSelected.Clone()
errorBand.SetFillColor(1)
errorBand.SetFillStyle(errorStyle)
HistSelected_1.SetFillStyle(errorStyle)
Data_hist.Draw("e1")
HistSelected.Draw("histsame")
HistSelected_1.Draw("e1same")

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
h3=HistSelected_1.Clone()
hwoE=HistSelected_1.Clone()
h1.SetMarkerStyle(20)
h1.SetMarkerSize(1.5)
#h1.Sumw2()
#h3.Sumw2()
h1.SetStats(0)
h1.Divide(hwoE)
h3.Divide(hwoE)
h3.SetMarkerSize(0.1)
h1.GetXaxis().SetTitle(xaxis_label)
h1.SetMarkerColor(1)
h1.SetLineColor(1)
h1.SetTitle("")
h1.GetXaxis().SetLabelSize(0.1)
h1.GetYaxis().SetTitle("ratio")
h1.GetYaxis().SetLabelSize(0.08)
h1.SetMaximum(2.5)
h1.SetMinimum(0.0)
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

c.Modified()
c.SaveAs("plots/plot_"+histoname+"_"+channel_+"_"+bkg_selected+"_nonorm.png")

