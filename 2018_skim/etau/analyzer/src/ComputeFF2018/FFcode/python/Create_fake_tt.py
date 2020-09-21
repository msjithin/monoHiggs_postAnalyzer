import ROOT
import re
from array import array
import argparse

def Create_fake_tt(directory):

   fileData=ROOT.TFile(directory+"/Data.root","r") #Data histogram
   fileVV=ROOT.TFile(directory+"/VV.root","r") #MC histogram with real leptons
   fileDY=ROOT.TFile(directory+"/DY.root","r") #MC histogram with real leptons
   fileTT=ROOT.TFile(directory+"/TT.root","r") #MC histogram with real leptons
   #fileW=ROOT.TFile(directory+"/W.root","r") #MC histogram with real leptons
   fileDataSub=ROOT.TFile(directory+"/Fake.root","recreate") #Data without real leptons
   
   dir_anti=["tt_mvis_anti","tt_met_anti","tt_tau1pt_anti","tt_tau2pt_anti","tt_pth_anti","tt_tau1eta_anti","tt_j1pt_anti"]
   dir_iso=["tt_mvis_iso","tt_met_iso","tt_tau1pt_iso","tt_tau2pt_iso","tt_pth_iso","tt_tau1eta_iso","tt_j1pt_iso"]
   ncat=7
   
   for i in range (0,ncat):
      Data=fileData.Get(dir_anti[i]).Get("data_obs")
      #Data.Add(fileVV.Get(dir_anti[i]).Get("VV"),-1.0)
      Data.Add(fileTT.Get(dir_anti[i]).Get("TT"),-1.0)
      DY=fileDY.Get(dir_anti[i]).Get("DY").Clone()
      Data.Add(DY,-1.0)
      Data.SetName("jetFakes")
      for k in range(1,Data.GetSize()-1):
        if Data.GetBinContent(k)<0:
   	   Data.SetBinContent(k,0)
      fileDataSub.cd()
      mydir=fileDataSub.mkdir(dir_iso[i])
      mydir.cd()
      Data.Write()
   
   
