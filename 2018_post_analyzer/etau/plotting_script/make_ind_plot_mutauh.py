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


parser = argparse.ArgumentParser()
parser.add_argument("-name",
                    help="name of hist to be plotted")

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
inputFile_=args.inputFile
category_=args.category
xaxis_label=args.xaxis
channel_=args.channel

#print "Item / histogram name = ", histoname
#print "histo 1               = ", histoname+"_"+category_
#print "histo 2               = ", histoname+"_"+category_

if (args.logYaxis):
  yaxisLog = 1
else :
  yaxisLog = 0

if (args.noLegend):
  noLegend = 1
elif (histoname =="muPt_Hpt_2D" or histoname =="muPt_Hpt_2D_highPt"):
  noLegend = 1
else :
  noLegend = 0
  
if (args.piAxis):
  piAxis = 1
else :
  piAxis = 0


def add_lumi():
    lowX=0.50
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.50, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   32 )#12
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.035)
    lumi.SetTextFont (   42 )
    if channel_=="combined":
      lumi.AddText("4 channels combined 2018, 100 fb^{-1} (13 TeV)")
    if channel_=="mutau":
      lumi.AddText("#mu#tau_{h} 2018, 100 fb^{-1} (13 TeV)")
    if channel_=="etau":
      lumi.AddText("e#tau_{h} 2018, 100 fb^{-1} (13 TeV)")
    if channel_=="tautau":
      lumi.AddText("#tau_{h}#tau_{h} 2018, 100 fb^{-1} (13 TeV)")
    if channel_=="emu":
      lumi.AddText("e#mu 2018, 100 fb^{-1} (13 TeV)")
    return lumi

def add_CMS():
    lowX=0.21
    lowY=0.70
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.08)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
    return lumi

def add_text_des(desctription, cat_select):
    lowX=0.50
    lowY=0.60
    lumi  = ROOT.TPaveText(lowX, lowY, lowX+0.40, lowY+0.2, "NDC")
    lumi.SetTextFont(42)
    lumi.SetTextSize(0.035)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    pt_select=''
    if (cat_select=='1') : 
      pt_select='0 <'+desctription+'<100'
    if (cat_select=='2') :
      pt_select='100 <'+desctription+'<200'
    if (cat_select=='3') :
      pt_select='200 <'+desctription+'<300'
    if (cat_select=='4') :
      pt_select='300 <'+desctription+'<400'
    if (cat_select=='5') :
      pt_select='400 <'+desctription+'<500'
    if (cat_select=='6') :
      pt_select='500 <'+desctription+'<600'
    if (cat_select=='7') :
      pt_select='600 <'+desctription+'<1000'
    lumi.AddText(pt_select)
    return lumi

def add_Preliminary():
    lowX=0.21
    lowY=0.63
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextFont (   40 )
    lumi.SetTextSize(0.06)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("Preliminary")
    return lumi

def make_legend():
        output = ROOT.TLegend(0.5, 0.75, 0.65, 0.84, "", "brNDC")
        #output = ROOT.TLegend(0.2, 0.1, 0.47, 0.65, "", "brNDC")
        output.SetLineWidth(1)
        output.SetLineStyle(1)
        output.SetFillStyle(0)
        output.SetBorderSize(1)
        output.SetTextFont(42)
        return output

#ROOT.gStyle.SetFrameLineWidth(1)
#ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetOptStat()
if histoname =="cutflow":
  ROOT.gStyle.SetOptStat(0)
if histoname =="cutflow_n" or histoname =="cutflow_Htt":
  ROOT.gStyle.SetOptStat(0)
if histoname =="muPt_Hpt_2D" or histoname =="muPt_Hpt_2D_highPt":
  ROOT.gStyle.SetOptStat(0)

c=ROOT.TCanvas("canvas","",0,0,1200,1200)
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

#OutFile_1=ROOT.TFile("ggh_aw_mutau.root","r")

#adapt=ROOT.gROOT.GetColor(12)
#new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
#trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)

Histo_org=OutFile.Get("Higgs_pt_gen_org")
Histo_ref=OutFile.Get("Higgs_pt_gen")
Histo_gen=OutFile.Get(histoname+"_"+category_)
ngen=OutFile.Get("nEvents").GetBinContent(1)
luminosity = 100000.0
xs = 44.14*0.0627
weight=luminosity*xs/ngen

if histoname=='cutflow' or histoname=='cutflow_n' or histoname =="cutflow_Htt" or histoname=='muPt_Hpt_2D' or histoname=='muPt_Hpt_2D_highPt' :
  Histo_gen=OutFile.Get(histoname)
  #weight=1.0

print "weight = ", weight
Histo_gen.Scale(weight)
Histo_org.Scale(weight)
Histo_ref.Scale(weight)
line_y_value=Histo_ref.Integral()
line_y_value_org=Histo_org.Integral()
myline=ROOT.TLine(0,line_y_value,11,line_y_value )
myline_org=ROOT.TLine(0,line_y_value_org,11,line_y_value_org )
myline.SetLineColor(4)
myline.SetLineWidth(4)
myline_org.SetLineColor(2)
myline_org.SetLineWidth(2)
#Histo_rec=OutFile.Get(histoname+"_"+category_)
#Histo_org.SetLineColor(1)
Histo_gen.SetMarkerStyle(20)
Histo_gen.SetMarkerSize(1.0)
Histo_gen.SetMarkerColor(4)
Histo_gen.SetFillColor(4)
Histo_gen.SetLineColor(4)
if histoname=='cutflow' or histoname=='cutflow_n' or histoname =="cutflow_Htt":
  Histo_gen.SetLineColor(46)
  Histo_gen.SetFillColor(46)
  Histo_gen.SetMarkerColor(46)
#Histo_rec.SetLineColor(2)
Histo_gen.SetTitle("")
if piAxis == 1:
  Histo_gen.GetXaxis().SetNdivisions(-405)
elif histoname=='Events_level_':
  Histo_gen.GetXaxis().SetNdivisions(-115)
elif (histoname=='mu_Charge' or histoname=='tau_Charge'):
  Histo_gen.GetXaxis().SetNdivisions(-104)
elif (histoname=='muCharge' or histoname=='tauCharge'):
  Histo_gen.GetXaxis().SetNdivisions(-108)
else :
  Histo_gen.GetXaxis().SetNdivisions(-405)

pad1 = ROOT.TPad("pad1","pad1",0,0.0,1,1)
pad1.Draw()
pad1.cd()
pad1.SetFillColor(0)
pad1.SetBorderMode(1)
pad1.SetBorderSize(2)
pad1.SetTickx(1)
pad1.SetTicky(1)
pad1.SetGridx()
pad1.SetLeftMargin(0.18)
pad1.SetRightMargin(0.05)
if histoname =="muPt_Hpt_2D" or histoname =="muPt_Hpt_2D_highPt":
  pad1.SetRightMargin(0.18)
pad1.SetTopMargin(0.122)
if histoname== "cutflow" or histoname=='cutflow_n' or histoname =="cutflow_Htt":
  pad1.SetBottomMargin(0.25)
else:
  pad1.SetBottomMargin(0.122)
pad1.SetFrameFillStyle(0)
pad1.SetFrameLineStyle(0)
pad1.SetFrameLineWidth(3)
pad1.SetFrameBorderMode(1)
pad1.SetFrameBorderSize(2)
if yaxisLog == 1 :
  pad1.SetLogy()  

Histo_gen.GetXaxis().SetLabelSize(0.04)
Histo_gen.GetYaxis().SetLabelFont(42)
Histo_gen.GetYaxis().SetLabelOffset(0.01)
Histo_gen.GetYaxis().SetLabelSize(0.04)
Histo_gen.GetYaxis().SetTitleSize(0.04)
Histo_gen.GetYaxis().SetTitleOffset(2.0)
Histo_gen.SetTitle("")
Histo_gen.GetYaxis().SetTitle("Events")
if( histoname=='muPt_Hpt_2D' or histoname=='muPt_Hpt_2D_highPt'):
  Histo_gen.GetYaxis().SetTitle("Higgs pt [GeV]")
if( histoname=='Hpt_muPt_2D'):                                                                                 
  Histo_gen.GetYaxis().SetTitle("mu pt [GeV]")    
Histo_gen.GetXaxis().SetTitle(xaxis_label)

if yaxisLog == 1 :
  Histo_gen.SetMaximum(10*Histo_gen.GetMaximum())
  if histoname== "cutflow" or histoname=='cutflow_n' or histoname =="cutflow_Htt": 
    Histo_gen.SetMinimum(100)
  else :
    Histo_gen.SetMinimum(0.01)

else :
  Histo_gen.SetMaximum(1.35*Histo_gen.GetMaximum())
  Histo_gen.SetMinimum(0.0)
#stack.SetLineColor(9)
#Histo_gen.Draw()
#Histo_rec.Draw("e1same")
#Histo_org.Draw("e1same")

#errorBand.Draw("e2same")
#ZH125_mc.Draw("histsame")

#Histo_gen.Draw("e1same")
#errorBand.Draw("e2same")
if histoname=='cutflow':
  Histo_gen.GetXaxis().SetLabelOffset(0.01)
  Histo_gen.GetXaxis().SetBinLabel(11,"#tau fake ejection");
  Histo_gen.GetXaxis().SetBinLabel(10,"deltaR");
  Histo_gen.GetXaxis().SetBinLabel( 9,"bjetVeto");
  Histo_gen.GetXaxis().SetBinLabel( 8,"ThirdLepVeto");
  Histo_gen.GetXaxis().SetBinLabel( 7,"MuTau charge");
  Histo_gen.GetXaxis().SetBinLabel( 6,"TauDecayMode");
  Histo_gen.GetXaxis().SetBinLabel( 5,"TauIso");
  #Histo_gen.GetXaxis().SetBinLabel( 8,"MuonIso");
  #Histo_gen.GetXaxis().SetBinLabel( 7,"MuonId");
  #Histo_gen.GetXaxis().SetBinLabel( 6,"MuonD0");
  #Histo_gen.GetXaxis().SetBinLabel( 5,"MuonDz");
  Histo_gen.GetXaxis().SetBinLabel( 4,"Tau selection");
  Histo_gen.GetXaxis().SetBinLabel( 3,"Muon selection");
  Histo_gen.GetXaxis().SetBinLabel( 2,"Trigger");
  Histo_gen.GetXaxis().SetBinLabel( 1,"Initial");
  Histo_gen.GetXaxis().SetNdivisions(-115)
  Histo_gen.GetXaxis().LabelsOption("v")
  Histo_gen.GetXaxis().SetTitle(" cut flow")
  Histo_gen.GetXaxis().SetTitleOffset(4)
  #  Histo_gen.SetMaximum(1.2)#FIXME(1.2)
  #  Histo_gen.SetMinimum(0.8)#FIXME(0.8)


if histoname=='cutflow_n':
  Histo_gen.GetXaxis().SetLabelOffset(0.01)
  Histo_gen.GetXaxis().SetBinLabel( 9,"#tau fake ejection");
  Histo_gen.GetXaxis().SetBinLabel( 8,"deltaR");
  Histo_gen.GetXaxis().SetBinLabel( 7,"bjetVeto");
  Histo_gen.GetXaxis().SetBinLabel( 6,"ThirdLepVeto");
  Histo_gen.GetXaxis().SetBinLabel( 5,"opposite charge");
  Histo_gen.GetXaxis().SetBinLabel( 4,"Tau selection");
  Histo_gen.GetXaxis().SetBinLabel( 3,"Muon selection");
  Histo_gen.GetXaxis().SetBinLabel( 2,"Trigger");
  Histo_gen.GetXaxis().SetBinLabel( 1,"Initial");
  Histo_gen.GetXaxis().SetNdivisions(-115)
  Histo_gen.GetXaxis().LabelsOption("v")
  Histo_gen.GetXaxis().SetTitle(" cut flow")
  Histo_gen.GetXaxis().SetTitleOffset(4)
  #  Histo_gen.SetMaximum(1.2)#FIXME(1.2)
  #  Histo_gen.SetMinimum(0.8)#FIXME(0.8)


if histoname=='cutflow_Htt':
  Histo_gen.GetXaxis().SetLabelOffset(0.01)
  Histo_gen.GetXaxis().SetBinLabel( 11,"Tau id");
  Histo_gen.GetXaxis().SetBinLabel( 10,"opposite charge");
  Histo_gen.GetXaxis().SetBinLabel( 9,"gen jet removal");
  Histo_gen.GetXaxis().SetBinLabel( 8,"Tau pt cut");
  Histo_gen.GetXaxis().SetBinLabel( 7,"btag veto");
  Histo_gen.GetXaxis().SetBinLabel( 6,"Lepton separation");
  Histo_gen.GetXaxis().SetBinLabel( 5,"Anti lepton disc.");
  Histo_gen.GetXaxis().SetBinLabel( 4,"Trigger");
  Histo_gen.GetXaxis().SetBinLabel( 3,"Met filters");
  Histo_gen.GetXaxis().SetBinLabel( 2,"Eta cuts");
  Histo_gen.GetXaxis().SetBinLabel( 1,"Skimming");
  Histo_gen.GetXaxis().SetNdivisions(-115)
  Histo_gen.GetXaxis().LabelsOption("v")
  Histo_gen.GetXaxis().SetTitle(" cut flow")
  Histo_gen.GetXaxis().SetTitleOffset(4)
  Histo_gen.SetMaximum(Histo_gen.GetMaximum()*1.2)#FIXME(1.2)
  Histo_gen.SetMinimum(Histo_gen.GetMinimum()*0.5)#FIXME(0.8)

if (histoname=='muPt_Hpt_2D' or histoname=='muPt_Hpt_2D_highPt'):
  Histo_gen.Draw("COLZ")
else:
  Histo_gen.Draw()
#if histoname=='cutflow':
#  myline.Draw()
#  myline_org.Draw()
legende=make_legend()
if histoname=="cutflow" or histoname=='cutflow_n' or histoname =="cutflow_Htt":
  legende.AddEntry(Histo_gen,"reco","f")
elif (category_ != "gen"):
  legende.AddEntry(Histo_gen,"reco","elp")
else :
  legende.AddEntry(Histo_gen,"gen","elp")

#legende.AddEntry(Histo_rec,"reco","elp")

if noLegend == 0 :
  legende.Draw()

l1=add_lumi()
l1.Draw("same")
l2=add_CMS()
#l2.Draw("same")
l3=add_Preliminary()
#l3.Draw("same")
l4=add_text_des('Higgs pt', category_)
if histoname=="muPt_Hpt" or histoname=="muEta_Hpt" :
  l4.Draw("same")
#pad1.RedrawAxis()

c.cd()
#ROOT.gPad.RedrawAxis()
c.Modified() 
c.SaveAs("plot_"+histoname+"_"+category_+".png") 
