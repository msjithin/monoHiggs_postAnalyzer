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

now = datetime.datetime.now()

firstarg=argv[1]
secondarg=argv[2]
thirdarg = argv[3]
fourtharg = ""
fiftharg = ""
if len(argv) > 4 :
  fourtharg=argv[4]

print "length = ",len(argv) 
thirdarg = thirdarg + " " + fourtharg
if firstarg != 'Events_level_' :
  histoname = firstarg+secondarg
else :
  histoname = firstarg
print "Item / histogram name = ", histoname
print "thirdarg = ", thirdarg

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts


myargs = getopts(argv)
if '-logYaxis' in myargs:  # Example usage.
  yaxisLog = 1
else :
  yaxisLog = 0

if '-noLegend' in myargs:  # Example usage.
  noLegend = _mc
else :
  noLegend = 0

if (firstarg == "h_dPhi_" or firstarg =="Tau_phi_" or firstarg =="Muon_phi_" or firstarg =="dphi_muMet_" or firstarg =="dphi_tauMet"):
  piAxis = 1
else :
  piAxis = 0



def add_lumi():
    lowX=0.50
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.10, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.06)
    lumi.SetTextFont (   42 )
    lumi.AddText("#mu#tau_{h}   2017, 41.52 fb^{-1} (13 TeV)")
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

def add_Preliminary():
    lowX=0.21
    lowY=0.63
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(52)
    lumi.SetTextSize(0.06)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("Preliminary")
    return lumi

def make_legend():
        output = ROOT.TLegend(0.70, 0.55, 0.92, 0.84, "", "brNDC")
        #output = ROOT.TLegend(0.2, 0.1, 0.47, 0.65, "", "brNDC")
        output.SetLineWidth(1)
        output.SetLineStyle(1)
        output.SetFillStyle(0)
        output.SetBorderSize(1)
        output.SetTextFont(62)
        return output

ROOT.gStyle.SetFrameLineWidth(1)
ROOT.gStyle.SetLineWidth(1)
ROOT.gStyle.SetOptStat(0)

c=ROOT.TCanvas("canvas","",0,0,600,600)
c.cd()

OutFile=ROOT.TFile("f_mutau_initial.root","r")
OutFile_fbkg=ROOT.TFile("ff_method/f_mutau_fakebkg.root","r")

adapt=ROOT.gROOT.GetColor(12)
new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)

categories=[firstarg,firstarg] 
dirName=[firstarg+"/","qcd/","bkg/"+firstarg+"/"] 
ncat=1


if firstarg=="MET" :  # Example usage.
  dirName[0]="MET_0_/"
if firstarg=="Cutflow": 
  dirName[0]="Cutflow_/"

for i in range (0,ncat):
    HistoData=OutFile.Get(dirName[0]+"data_obs_"+histoname)
    
    #    ZTT=OutFile.Get(dirName[0]+"DY_LO_"+histoname)
    VV=OutFile.Get(dirName[0]+"VV_"+histoname)

    ZTT=OutFile.Get(dirName[0]+"ZTTinc_"+histoname)
    ZTT1=OutFile.Get(dirName[0]+"ZTT1jet_"+histoname)
    ZTT2=OutFile.Get(dirName[0]+"ZTT2jet_"+histoname)
    ZTT3=OutFile.Get(dirName[0]+"ZTT3jet_"+histoname)
    ZTT4=OutFile.Get(dirName[0]+"ZTT4jet_"+histoname)
    ZTT.Add(ZTT1)
    ZTT.Add(ZTT2)
    ZTT.Add(ZTT3)
    ZTT.Add(ZTT4)

    EWKZ_0=OutFile.Get(dirName[0]+"EWKZ2Jets_ZToLL_"+histoname)
    EWKZ=OutFile.Get(dirName[0]+"EWKZ2Jets_ZToNuNu_"+histoname)
    EWKZ.Add(EWKZ_0)
    EWKW_0=OutFile.Get(dirName[0]+"EWKWMinus2Jets_"+histoname)
    EWKW=OutFile.Get(dirName[0]+"EWKWPlus2Jets_"+histoname)
    EWKW.Add(EWKW_0)

    ZTT.Add(EWKZ)

    #TTJ=OutFile.Get(dirName[0]+"TTJets_"+histoname)
    TTJ=OutFile.Get(dirName[0]+"TTTo2L2Nu_"+histoname)
    TTJ_1=OutFile.Get(dirName[0]+"TTToHadronic_"+histoname)
    TTJ_2=OutFile.Get(dirName[0]+"TTToSemiLeptonic_"+histoname)
    TTJ.Add(TTJ_1)
    TTJ.Add(TTJ_2)

    ggH125=OutFile.Get(dirName[0]+"GluGluHToTauTau_"+histoname)
    qqH125=OutFile.Get(dirName[0]+"VBFHToTauTau_"+histoname)
    qqH125_1=OutFile.Get(dirName[0]+"VBFHToWWTo2L2Nu_"+histoname)
    ZHToTauTau=OutFile.Get(dirName[0]+"ZHToTauTau_"+histoname)
    WpH125=OutFile.Get(dirName[0]+"WPlusHToTauTau_"+histoname)
    WmH125=OutFile.Get(dirName[0]+"WMinusHToTauTau_"+histoname)
    ggH125.Add(ZHToTauTau)
    ggH125.Add(qqH125)
    ggH125.Add(qqH125_1)
    ggH125.Add(WpH125)
    ggH125.Add(WmH125)
    
    sT_t_antitop=OutFile.Get(dirName[0]+"ST_t-channel_antitop_"+histoname)
    sT_t_top=OutFile.Get(dirName[0]+"ST_t-channel_top_"+histoname)
    sT_tW_antitop=OutFile.Get(dirName[0]+"ST_tW_antitop_"+histoname)
    sT_tW_top=OutFile.Get(dirName[0]+"ST_tW_top_"+histoname)
    sT_t_antitop.Add(sT_t_top)
    sT_t_antitop.Add(sT_tW_antitop)
    sT_t_antitop.Add(sT_tW_top)

    VVV=OutFile.Get(dirName[0]+"VVV_"+histoname)
    ZJetsToNuNu=OutFile.Get(dirName[0]+"ZJetsToNuNu_"+histoname)
    WWJJ=OutFile.Get(dirName[0]+"WpWpJJ_EWK_QCD_"+histoname)
 #   ggWW=OutFile.Get(dirName[0]+"GluGluWWTo2L2Nu_"+histoname)
    VV.Add(WWJJ)
    VV.Add(VVV)
    VV.Add(sT_t_antitop)
  #  VV.Add(ggWW)

#    WJets=OutFile.Get(dirName[0]+"WJetsToLNu_"+histoname)
#    WJets_100To200=OutFile.Get(dirName[0]+"WJetsToLNu_HT100To200_"+histoname)
#    WJets_200To400=OutFile.Get(dirName[0]+"WJetsToLNu_HT200To400_"+histoname)
#    WJets_400To600=OutFile.Get(dirName[0]+"WJetsToLNu_HT400To600_"+histoname)
#    WJets_600To800=OutFile.Get(dirName[0]+"WJetsToLNu_HT600To800_"+histoname)
#    WJets_800To1200=OutFile.Get(dirName[0]+"WJetsToLNu_HT800To1200_"+histoname)
#    WJets_1200To2500=OutFile.Get(dirName[0]+"WJetsToLNu_HT1200To2500_"+histoname)
#    WJets_2500ToInf=OutFile.Get(dirName[0]+"WJetsToLNu_HT2500ToInf_"+histoname)
#    nbins_ =  WJets.GetNbinsX()
#    print "WJets inclusive events   ", WJets.Integral(0,nbins_+1)
#    print "WJets 100-200 events   ", WJets_100To200.Integral(0,nbins_+1)
#    print "WJets 200-400 events   ", WJets_200To400.Integral(0,nbins_+1)
#    print "WJets 400-600 events   ", WJets_400To600.Integral(0,nbins_+1)
#    print "WJets 600-800 events   ", WJets_600To800.Integral(0,nbins_+1)
#    print "WJets 800-1200 events   ", WJets_800To1200.Integral(0,nbins_+1)
#    print "WJets 1200-2500 events   ", WJets_1200To2500.Integral(0,nbins_+1)
#    print "WJets 2500-Inf events   ", WJets_2500ToInf.Integral(0,nbins_+1)

#    WJets.Add(WJets_100To200)
#    WJets.Add(WJets_200To400)
#    WJets.Add(WJets_400To600)
 #   WJets.Add(WJets_600To800)
 #   WJets.Add(WJets_800To1200)
 #   WJets.Add(WJets_1200To2500)
 #   WJets.Add(WJets_2500ToInf)
 #   WJets.Add(EWKW)

    F_bkg=OutFile_fbkg.Get(dirName[0]+"data_obs_"+histoname)
    ZTT_fbkg=OutFile_fbkg.Get(dirName[0]+"ZTTinc_"+histoname)
    VV_fbkg=OutFile_fbkg.Get(dirName[0]+"VV_"+histoname)
    ZTT1_fbkg=OutFile_fbkg.Get(dirName[0]+"ZTT1jet_"+histoname)
    ZTT2_fbkg=OutFile_fbkg.Get(dirName[0]+"ZTT2jet_"+histoname)
    ZTT3_fbkg=OutFile_fbkg.Get(dirName[0]+"ZTT3jet_"+histoname)
    ZTT4_fbkg=OutFile_fbkg.Get(dirName[0]+"ZTT4jet_"+histoname)
    EWKZ_0_fbkg=OutFile_fbkg.Get(dirName[0]+"EWKZ2Jets_ZToLL_"+histoname)
    EWKZ_fbkg=OutFile_fbkg.Get(dirName[0]+"EWKZ2Jets_ZToNuNu_"+histoname)
    EWKW_0_fbkg=OutFile_fbkg.Get(dirName[0]+"EWKWMinus2Jets_"+histoname)
    EWKW_fbkg=OutFile_fbkg.Get(dirName[0]+"EWKWPlus2Jets_"+histoname)
    
    TTJ_fbkg=OutFile_fbkg.Get(dirName[0]+"TTTo2L2Nu_"+histoname)
    TTJ_fbkg1=OutFile_fbkg.Get(dirName[0]+"TTToHadronic_"+histoname)
    TTJ_fbkg2=OutFile_fbkg.Get(dirName[0]+"TTToSemiLeptonic_"+histoname)
    
    ggH125_fbkg=OutFile_fbkg.Get(dirName[0]+"GluGluHToTauTau_"+histoname)
    qqH125_fbkg=OutFile_fbkg.Get(dirName[0]+"VBFHToTauTau_"+histoname)
    ZHToTauTau_fbkg=OutFile_fbkg.Get(dirName[0]+"ZHToTauTau_"+histoname)
    WpH125_fbkg=OutFile_fbkg.Get(dirName[0]+"WPlusHToTauTau_"+histoname)
    WmH125_fbkg=OutFile_fbkg.Get(dirName[0]+"WMinusHToTauTau_"+histoname)
    sT_t_antitop_fbkg=OutFile_fbkg.Get(dirName[0]+"ST_t-channel_antitop_"+histoname)
    sT_t_top_fbkg=OutFile_fbkg.Get(dirName[0]+"ST_t-channel_top_"+histoname)
    sT_tW_antitop_fbkg=OutFile_fbkg.Get(dirName[0]+"ST_tW_antitop_"+histoname)
    sT_tW_top_fbkg=OutFile_fbkg.Get(dirName[0]+"ST_tW_top_"+histoname)
    VVV_fbkg=OutFile_fbkg.Get(dirName[0]+"VVV_"+histoname)
    ZJetsToNuNu_fbkg=OutFile_fbkg.Get(dirName[0]+"ZJetsToNuNu_"+histoname)
    WWJJ_fbkg=OutFile_fbkg.Get(dirName[0]+"WpWpJJ_EWK_QCD_"+histoname)
    
    F_bkg.Add(ZTT_fbkg, -1)
    F_bkg.Add(VV_fbkg, -1)
    F_bkg.Add(ZTT1_fbkg, -1)
    F_bkg.Add(ZTT2_fbkg, -1)
    F_bkg.Add(ZTT3_fbkg, -1)
    F_bkg.Add(ZTT4_fbkg, -1)
    F_bkg.Add(EWKZ_0_fbkg, -1)
    F_bkg.Add(EWKZ_fbkg, -1)
    F_bkg.Add(EWKW_0_fbkg, -1)
    F_bkg.Add(EWKW_fbkg, -1)
    F_bkg.Add(TTJ_fbkg, -1)
    F_bkg.Add(TTJ_fbkg1, -1)   
    F_bkg.Add(TTJ_fbkg2, -1)   
    F_bkg.Add(ggH125_fbkg, -1)
    F_bkg.Add(qqH125_fbkg, -1)
    F_bkg.Add(ZHToTauTau_fbkg, -1)
    F_bkg.Add(WpH125_fbkg, -1)
    F_bkg.Add(WmH125_fbkg, -1)
    F_bkg.Add(sT_t_antitop_fbkg, -1)
    F_bkg.Add(sT_t_top_fbkg, -1)
    F_bkg.Add(sT_tW_antitop_fbkg, -1)
    F_bkg.Add(sT_tW_top_fbkg, -1)
    F_bkg.Add(VVV_fbkg, -1)
    F_bkg.Add(ZJetsToNuNu_fbkg, -1)
    F_bkg.Add(WWJJ_fbkg, -1)



    ZTT.SetFillColor(ROOT.TColor.GetColor("#ffc866"))
    VV.SetFillColor(ROOT.TColor.GetColor("#9600a0 "))
    TTJ.SetFillColor(ROOT.TColor.GetColor("#9899cf"))
    ggH125.SetFillColor(ROOT.TColor.GetColor("#0165ff"))
    F_bkg.SetFillColor(ROOT.TColor.GetColor("#ce6468"))
#    WWW.SetFillColor(ROOT.TColor.GetColor("#99ccff"))
    ZJetsToNuNu.SetFillColor(ROOT.TColor.GetColor("#fe9a9c"))


    #EWKZ.SetLineColor(1)
    ZTT.SetLineColor(1)
    VV.SetLineColor(1)
    TTJ.SetLineColor(1)
    ggH125.SetLineColor(1)
    F_bkg.SetLineColor(1)
 #   WWW.SetLineColor(1)
    ZJetsToNuNu.SetLineColor(1)
    #QCD_mc.SetLineColor(1)
   # ZH125_mc.SetLineColor(4)
  #  ZH125_mc.SetLineWidth(5)

    stack=ROOT.THStack("stack","stack")
    stack.Add(ZJetsToNuNu)
  #  stack.Add(WWW)
    stack.Add(ggH125)
    stack.Add(VV)
    stack.Add(F_bkg)
    stack.Add(TTJ)
    stack.Add(ZTT)
   
    if histoname=='Events_level_':
      print "Initial "
      print "DY events = ", ZTT.GetBinContent(1)
      print "VV events = ", VV.GetBinContent(1)
      print "TT events = ", TTJ.GetBinContent(1)
      print "F_bkg events = ", F_bkg.GetBinContent(1)
      print "Data events = ",HistoData.GetBinContent(1)
      
      print "After b-jet veto "
      print "DY events = ", ZTT.GetBinContent(9)
      print "VV events = ", VV.GetBinContent(9)
      print "TT events = ", TTJ.GetBinContent(9)
      print "F_bkg events = ", F_bkg.GetBinContent(9)
      print "Data events = ",HistoData.GetBinContent(9)
      text_file = open("Output.txt", "a+")
      text_file.write("\n")
      text_file.write("*************  From cut flow *************" )
      text_file.write("\n")
      text_file.write("Initial" )
      text_file.write("\n")
      text_file.write("DY events = "+ str(ZTT.GetBinContent(1))) 
      text_file.write("\n")
      text_file.write("VV events = "+ str(VV.GetBinContent(1) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(1) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(1) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(1)))
      text_file.write("\n")
      text_file.write("********************************")
      text_file.write("\n")
      text_file.write(" passed MET filter" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(2) ))
      text_file.write("\n")
      text_file.write("VV events = "+str( VV.GetBinContent(2) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(2) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(2) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(2)))
      text_file.write("\n")
      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("ele trigger" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(3) ))
      text_file.write("\n")
      text_file.write("VV events = "+str( VV.GetBinContent(3) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(3) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(3) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(3)))
      text_file.write("\n")



      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("good ele" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(4) ))
      text_file.write("\n")
      text_file.write("VV events = "+str( VV.GetBinContent(4) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(4) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(4) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(4)))
      text_file.write("\n")

      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("good tau" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(5) ))
      text_file.write("\n")
      text_file.write("VV events = "+str( VV.GetBinContent(5) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(5) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(5) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(5)))
      text_file.write("\n")
      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("opp charge" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(6) ))
      text_file.write("\n")
      text_file.write("VV events = "+str( VV.GetBinContent(6) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(6) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(6) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(6)))
      text_file.write("\n")

      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("delta R" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(7) ))
      text_file.write("\n")
      text_file.write("VV events = "+str( VV.GetBinContent(7) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(7) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(7) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(7)))
      text_file.write("\n")
      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("3rd lepton veto" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(8) ))
      text_file.write("\n")
      text_file.write("VV events = "+str( VV.GetBinContent(8) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(8) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(8) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(8)))
      text_file.write("\n")

      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("After b-jet veto" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(9) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(9) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(9) ))
      text_file.write("\n")
      text_file.write(" VV events = "+str( VV.GetBinContent(9) ))
      text_file.write("\n")
      text_file.write(" ggH125 events = "+str( ggH125.GetBinContent(9) ))
      text_file.write("\n")
      text_file.write(" ZJetsToNuNu events = "+str( ZJetsToNuNu.GetBinContent(9) ))
      text_file.write("\n")
      text_file.write(" Total bkg = "+str( ZTT.GetBinContent(9)+TTJ.GetBinContent(9)+F_bkg.GetBinContent(9)+VV.GetBinContent(9)+ggH125.GetBinContent(9)+ZJetsToNuNu.GetBinContent(9) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(9)))
      text_file.write("\n")
      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("After Higgs pt cut" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(10) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(10) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(10) ))
      text_file.write("\n")
      text_file.write(" VV events = "+str( VV.GetBinContent(10) ))
      text_file.write("\n")
      text_file.write(" ggH125 events = "+str( ggH125.GetBinContent(10) ))
      text_file.write("\n")
      text_file.write(" ZJetsToNuNu events = "+str( ZJetsToNuNu.GetBinContent(10) ))
      text_file.write("\n")
      text_file.write(" Total bkg = "+str( ZTT.GetBinContent(10)+TTJ.GetBinContent(10)+F_bkg.GetBinContent(10)+VV.GetBinContent(10)+ggH125.GetBinContent(10)+ZJetsToNuNu.GetBinContent(10) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(10)))

      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("After visible mass cut" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(11) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(11) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(11) ))
      text_file.write("\n")
      text_file.write(" VV events = "+str( VV.GetBinContent(11) ))
      text_file.write("\n")
      text_file.write(" ggH125 events = "+str( ggH125.GetBinContent(11) ))
      text_file.write("\n")
      text_file.write(" ZJetsToNuNu events = "+str( ZJetsToNuNu.GetBinContent(11) ))
      text_file.write("\n")
      text_file.write(" Total bkg = "+str( ZTT.GetBinContent(11)+TTJ.GetBinContent(11)+F_bkg.GetBinContent(11)+VV.GetBinContent(11)+ggH125.GetBinContent(11)+ZJetsToNuNu.GetBinContent(11) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(11)))
      text_file.write("\n")
      text_file.write("********************************")
      text_file.write("\n")
      text_file.write("After met cut" )
      text_file.write("\n")
      text_file.write("DY events = "+str(ZTT.GetBinContent(12) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.GetBinContent(12) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.GetBinContent(12) ))
      text_file.write("\n")
      text_file.write(" VV events = "+str( VV.GetBinContent(12) ))
      text_file.write("\n")
      text_file.write(" ggH125 events = "+str( ggH125.GetBinContent(12) ))
      text_file.write("\n")
      text_file.write(" ZJetsToNuNu events = "+str( ZJetsToNuNu.GetBinContent(12) ))
      text_file.write("\n")
      text_file.write("\n")
      text_file.write(" Total bkg = "+str( ZTT.GetBinContent(12)+TTJ.GetBinContent(12)+F_bkg.GetBinContent(12)+VV.GetBinContent(12)+ggH125.GetBinContent(12)+ZJetsToNuNu.GetBinContent(12) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.GetBinContent(12)))
      text_file.write("\n")

      text_file.close()
      
      
      
#      text_file.close()
      
    if histoname=='Electron_Pt_3' or histoname=='h_dPhi_3' or histoname=='Electron_phi_3'or histoname=='Tau_phi_3' or histoname=='Electron_eta_3'or histoname=='Tau_eta_3':
        
      nbins_ZTT =  ZTT.GetNbinsX()
      nbins_VV =  VV.GetNbinsX()
      nbins_TT =  TTJ.GetNbinsX()
      nbins_F_bkg =  F_bkg.GetNbinsX()
      nbins_Data =  HistoData.GetNbinsX()
      
      print "sameple ", histoname
      print "After b-jet veto "
      print "DY events = ", ZTT.Integral(0,nbins_ZTT+1) 
      print "VV events = ", VV.Integral(0,nbins_VV+1) 
      print "TT events = ", TTJ.Integral(0, nbins_TT+1)
      print "F_bkg events = ", F_bkg.Integral(0,nbins_F_bkg+1) 
      print "Data events = ",HistoData.Integral(0,nbins_Data+1) 
      
      text_file = open("Output.txt", "a")
      text_file.write("\n")
      text_file.write("**************  histogram name = "+ histoname + "************" )
      text_file.write("\n")
      text_file.write("After b-jet veto" )
      text_file.write("\n")
      text_file.write("DY events =  "+ str( ZTT.Integral(0,nbins_ZTT+1) ))
      text_file.write("\n")
      text_file.write("VV events = "+ str( VV.Integral(0,nbins_VV+1) ))
      text_file.write("\n")
      text_file.write("TT events = "+str(TTJ.Integral(0,nbins_TT+1) ))
      text_file.write("\n")
      text_file.write("F_bkg events = "+str(F_bkg.Integral(0,nbins_F_bkg+1) ))
      text_file.write("\n")
      text_file.write("Data events = "+str( HistoData.Integral(0,nbins_Data+1) ))
      text_file.write("\n")
      text_file.write("\n")
      
      text_file.close()
    



    HistoData.GetXaxis().SetTitle("")
    HistoData.GetXaxis().SetTitleSize(0)


    if piAxis == 1:
      HistoData.GetXaxis().SetNdivisions(-405)
    elif histoname=='Events_level_':
      HistoData.GetXaxis().SetNdivisions(-115)
    else :
      HistoData.GetXaxis().SetNdivisions(-405)

    HistoData.GetXaxis().SetLabelSize(0)
    HistoData.GetYaxis().SetLabelFont(42)
    HistoData.GetYaxis().SetLabelOffset(0.01)
    HistoData.GetYaxis().SetLabelSize(0.05)
    HistoData.GetYaxis().SetTitleSize(0.075)
    HistoData.GetYaxis().SetTitleOffset(1.04)
    HistoData.SetTitle("")
    HistoData.GetYaxis().SetTitle("Events")

    errorBand=ZTT.Clone()
#    errorBand=OutFile.Get(dirName[0]+"DY_LO_"+histoname).Clone()
#    errorBand.Add(ZTT1)
#    errorBand.Add(ZTT2)
#    errorBand.Add(ZTT3)
#    errorBand.Add(ZTT4)
#    errorBand.Add(EWKZ)
    errorBand.Add(VV)
    errorBand.Add(TTJ)
    errorBand.Add(ggH125)
    errorBand.Add(F_bkg)
    #errorBand.Add(EWKZ)
   # errorBand.Add(WWW)
    errorBand.Add(ZJetsToNuNu)

    g=1
    while (g<errorBand.GetNbinsX()+1 and histoname !="Events_level_" and histoname !="Cutflow_"):
      ZTTerror = ZTT.GetBinContent(g)*0.05
      TTJerror = TTJ.GetBinContent(g)*0.06
      F_bkgerror = F_bkg.GetBinContent(g)*0.3
      VVerror = VV.GetBinContent(g)*0.05
      ggH125error = ggH125.GetBinContent(g)*0.05
      ZJetsToNuNuerror = ZJetsToNuNu.GetBinContent(g)*0.1
      total_err = sqrt((ZTTerror*ZTTerror) + (TTJerror*TTJerror) + (F_bkgerror*F_bkgerror) + (VVerror*VVerror) + (ggH125error*ggH125error) + (ZJetsToNuNuerror*ZJetsToNuNuerror) )
      errorBand.SetBinError(g, total_err)
      g = g+1


    errorBand.SetMarkerSize(0)
    errorBand.SetFillColor(1)
    errorBand.SetFillStyle(3001)
    errorBand.SetLineWidth(1)

pad1 = ROOT.TPad("pad1","pad1",0,0.35,1,1)
pad1.Draw()
pad1.cd()
pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetBorderSize(1)
pad1.SetTickx(1)
pad1.SetTicky(1)
pad1.SetGridx()
pad1.SetLeftMargin(0.18)
pad1.SetRightMargin(0.05)
pad1.SetTopMargin(0.122)
pad1.SetBottomMargin(0.026)
pad1.SetFrameFillStyle(0)
pad1.SetFrameLineStyle(0)
pad1.SetFrameLineWidth(1)
pad1.SetFrameBorderMode(0)
pad1.SetFrameBorderSize(1)
if yaxisLog == 1 :
  pad1.SetLogy()  

#for k in range(1,Data.GetSize()-1):
    #s=ZH125.GetBinContent(k)
    ##b=VV.GetBinContent(k)+Fake.GetBinContent(k)
    #b=VV.GetBinContent(k)+ggH125.GetBinContent(k)+W.GetBinContent(k)
    #if (b<0):
	#b=0.000001
    #if (s/(0.00001+0.05*s+b)**0.5 > 0.8):
	#Data.SetBinContent(k,-1)
	#Data.SetBinError(k,-1)
HistoData.SetMarkerStyle(20)

HistoData.SetMarkerSize(0.75)
HistoData.GetXaxis().SetTitle(thirdarg)
HistoData.GetYaxis().SetTitle("Events")
if yaxisLog == 1 :
  HistoData.SetMaximum(1000*max(HistoData.GetMaximum(),stack.GetMaximum()))
  HistoData.SetMinimum(0.1)
else :
  HistoData.SetMaximum(1.35*max(HistoData.GetMaximum(),stack.GetMaximum()))
  HistoData.SetMinimum(0.0)
#stack.SetLineColor(9)
HistoData.Draw("e1")
stack.Draw("histsame")

errorBand.Draw("e2same")
#ZH125_mc.Draw("histsame")

HistoData.Draw("e1same")
#errorBand.Draw("e2same")

legende=make_legend()
legende.AddEntry(HistoData,"Data","elp")
#legende.AddEntry(HistoData,"Data runB","elp")
#legende.AddEntry(HistoData,"Data runC","elp")
#legende.AddEntry(HistoData,"Data runD","elp")
#legende.AddEntry(HistoData,"Data runE","elp")
#legende.AddEntry(HistoData,"Data runF","elp")

#legende.AddEntry(HistoData,"fake background","f")
#legende.AddEntry(W_mc,"W","f")
legende.AddEntry(ZTT,"Drell-Yan","f")
legende.AddEntry(TTJ,"TT jets","f")
legende.AddEntry(VV,"VV, single top, VVV","f")
#legende.AddEntry(QCD_mc,"QCD","f")
legende.AddEntry(ggH125,"ggH, VBF, WH, ZH","f")
legende.AddEntry(F_bkg,"fake background","f")
#legende.AddEntry(EWKZ,"EWK","f")
#legende.AddEntry(WWW,"VVV","f")
legende.AddEntry(ZJetsToNuNu,"ZJetsToNuNu","f")

if noLegend == 0:
  legende.Draw()

l1=add_lumi()
l1.Draw("same")
l2=add_CMS()
l2.Draw("same")
l3=add_Preliminary()
l3.Draw("same")

pad1.RedrawAxis()

categ  = ROOT.TPaveText(0.21, 0.5+0.013, 0.43, 0.70+0.155, "NDC")
categ.SetBorderSize(   0 )
categ.SetFillStyle(    0 )
categ.SetTextAlign(   12 )
categ.SetTextSize ( 0.06 )
categ.SetTextColor(    1 )
categ.SetTextFont (   42 )
if i+1==1:       
  categ.AddText("OS")
elif i+1==2:       
  categ.AddText("SS")
#categ.Draw()

c.cd()
pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.35);
pad2.SetTopMargin(0.05);
pad2.SetBottomMargin(0.35);
pad2.SetLeftMargin(0.18);
pad2.SetRightMargin(0.05);
pad2.SetTickx(1)
pad2.SetTicky(1)
pad2.SetFrameLineWidth(1)
#pad2.SetGridx()
pad2.SetGridy()
pad2.Draw()
pad2.cd()
h1=HistoData.Clone()
h1.SetMaximum(2.0)#FIXME(1.5)
h1.SetMinimum(0.0)#FIXME(0.5)
h1.SetMarkerStyle(20)
h3=errorBand.Clone()
hwoE=errorBand.Clone()
for iii in range (1,hwoE.GetSize()-2):
  hwoE.SetBinError(iii,0)
h3.Sumw2()
h1.Sumw2()
h1.SetStats(0)
h1.Divide(hwoE)
h3.Divide(hwoE)
h1.GetXaxis().SetTitle(thirdarg)#(#vec{p_{T}}(#tau_{1})+#vec{p_{T}}(#tau_{2}))/(p_{T}(#tau_{1})+p_{T}(#tau_{2}))")#("m_{vis} (GeV)")#(#vec{p_{T}(#mu)}+#vec{p_{T}(#tau)})/(p_{T}(#mu)+p_{T}(#tau))")
#if (i+1==1 or i+1==2 or i+1==7 or i+1==8):
#	h1.GetXaxis().SetTitle("Electron p_{T} (GeV)")
#if (i+1==4 or i+1==10):
#     h1.GetXaxis().SetTitle("Muon p_{T} (GeV)")
#if (i+1==6 or i+1==12 or i+1==3 or i+1==5 or i+1==9 or i+1==11):
#     h1.GetXaxis().SetTitle("Tau p_{T} (GeV)")
h1.GetXaxis().SetLabelSize(0.1)
h1.GetYaxis().SetTitle("Obs./Exp.")
h1.GetYaxis().SetLabelSize(0.11)
h1.SetMaximum(1.5)#FIXME(1.5)
h1.SetMinimum(0.5)#FIXME(0.5)
if piAxis ==1 :
  if firstarg == "Tau_phi_" or firstarg == "Muon_phi_" :
    h1.GetXaxis().SetBinLabel(21,"#pi");
    h1.GetXaxis().SetBinLabel(16,"#frac{#pi}{2}");
    h1.GetXaxis().SetBinLabel(11,"0");
    h1.GetXaxis().SetBinLabel(6,"#frac{-#pi}{2}");
    h1.GetXaxis().SetBinLabel(1,"-#pi");
  else :
    h1.GetXaxis().SetBinLabel(20,"#pi");
    h1.GetXaxis().SetBinLabel(15,"#frac{3#pi}{4}");
    h1.GetXaxis().SetBinLabel(10,"#frac{#pi}{2}");
    h1.GetXaxis().SetBinLabel(5,"#frac{#pi}{4}");
    h1.GetXaxis().SetBinLabel(1,"0");

if piAxis == 1:
  h1.GetXaxis().SetLabelOffset(0.01)
  h1.GetXaxis().SetLabelSize(0.15)
  h1.GetXaxis().SetNdivisions(-405)

elif histoname=='Events_level_':
  
  h1.GetXaxis().SetBinLabel(12,"METcut");
  h1.GetXaxis().SetBinLabel(11,"VisibleMass");
  h1.GetXaxis().SetBinLabel(10,"HiggsPt");
  h1.GetXaxis().SetBinLabel(9,"Mt cut");
  h1.GetXaxis().SetBinLabel(8,"BjetVeto");
  h1.GetXaxis().SetBinLabel(7,"ThirdLepVeto");
  h1.GetXaxis().SetBinLabel(6,"DeltaR");
  h1.GetXaxis().SetBinLabel(5,"tau &opp. charge");
  #h1.GetXaxis().SetBinLabel(5,"GoodTau");
  h1.GetXaxis().SetBinLabel(4,"GoodMuon");
  h1.GetXaxis().SetBinLabel(3,"mu Trg");
  h1.GetXaxis().SetBinLabel(2,"MET Filters");
  h1.GetXaxis().SetBinLabel(1,"No cuts");
  h1.GetXaxis().SetNdivisions(-115)
  h1.GetXaxis().LabelsOption("v")
  h1.GetXaxis().SetTitle(" cut flow")
  h1.SetMaximum(1.2)#FIXME(1.2)
  h1.SetMinimum(0.8)#FIXME(0.8)

elif histoname=='Cutflow_':

  h1.GetXaxis().SetBinLabel(5,"BjetVeto");
  h1.GetXaxis().SetBinLabel(4,"ThirdLepVeto");
  h1.GetXaxis().SetBinLabel(3,"DeltaR");
  h1.GetXaxis().SetBinLabel(2,"tau & opp. charge");
  #h1.GetXaxis().SetBinLabel(2,"GoodTau");
  h1.GetXaxis().SetBinLabel(1,"GoodMu");
  h1.GetXaxis().SetNdivisions(-115)
  h1.GetXaxis().LabelsOption("v")
  h1.GetXaxis().SetTitle(" cut flow")
  h1.SetMaximum(1.2)#FIXME(1.2)
  h1.SetMinimum(0.8)#FIXME(0.8)

else :
  h1.GetXaxis().SetNdivisions(-405)

  
h1.GetYaxis().SetNdivisions(5)
h1.GetXaxis().SetTitleFont(42)
h1.GetYaxis().SetTitleFont(42)


h1.GetXaxis().SetTitleSize(0.15)
h1.GetYaxis().SetTitleSize(0.15)
h1.GetYaxis().SetTitleOffset(0.56)
h1.GetXaxis().SetTitleOffset(1.1)






h1.Draw("e0p")
h3.Draw("e2same")

c.cd()
#  pad1.Draw()

ROOT.gPad.RedrawAxis()

c.Modified()
#c.SaveAs("plots/"+firstarg+"_fr.pdf")
c.SaveAs("plots/"+firstarg+secondarg+".png")


