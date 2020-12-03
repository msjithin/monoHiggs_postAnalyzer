import ROOT
import re
from array import array
import argparse

def Subtract_prompt_tt(directory):

   fileData=ROOT.TFile(directory+"/Data.root","r") #Data histogram
   fileVV=ROOT.TFile(directory+"/VV.root","r") #MC histogram with real leptons
   fileDY=ROOT.TFile(directory+"/DY.root","r") #MC histogram with real leptons
   fileTT=ROOT.TFile(directory+"/TT.root","r") #MC histogram with real leptons
   #fileW=ROOT.TFile(directory+"/W.root","r") #MC histogram with real leptons
   fileDataSub=ROOT.TFile(directory+"/DataSub.root","recreate") #Data without real leptons
   
   dir_qcd=["tt_0jet_qcd_anti","tt_0jet_qcd_iso","tt_0SSloose_qcd_anti","tt_0SSloose_qcd_iso","tt_1jet_qcd_anti","tt_1jet_qcd_iso","tt_1SSloose_qcd_anti","tt_1SSloose_qcd_iso","tt_2jet_qcd_anti","tt_2jet_qcd_iso","tt_2SSloose_qcd_anti","tt_2SSloose_qcd_iso"]
   ncat_qcd=12
   
   if "corr1" in directory:
     dir_qcd=["tt_0jet_qcd_mvis_anti","tt_0jet_qcd_mvis_iso","tt_0SSloose_qcd_mvis_anti","tt_0SSloose_qcd_mvis_iso","tt_0jet_qcd_tau2pt_anti","tt_0jet_qcd_tau2pt_iso","tt_0SSloose_qcd_tau2pt_anti","tt_0SSloose_qcd_tau2pt_iso","tt_0jet_qcd_tau1pt_anti","tt_0jet_qcd_tau1pt_iso","tt_0SSloose_qcd_tau1pt_anti","tt_0SSloose_qcd_tau1pt_iso","tt_0jet_qcd_met_anti","tt_0jet_qcd_met_iso","tt_0SSloose_qcd_met_anti","tt_0SSloose_qcd_met_iso"]
     ncat_qcd=16
   
   if "OSSSFF" in directory:
     ncat_qcd=2
   
   for i in range (0,ncat_qcd):
      Data=fileData.Get(dir_qcd[i]).Get("data_obs")
      Data.Add(fileVV.Get(dir_qcd[i]).Get("VVLT"),-1.0)
      Data.Add(fileVV.Get(dir_qcd[i]).Get("VVJ"),-1.0)
      Data.Add(fileTT.Get(dir_qcd[i]).Get("TTJ"),-1.0)
      Data.Add(fileTT.Get(dir_qcd[i]).Get("TTLT"),-1.0)
      Data.Add(fileDY.Get(dir_qcd[i]).Get("DYJ"),-1.0)
      DYLT=fileDY.Get(dir_qcd[i]).Get("DYLT").Clone()
      for j in range(1,DYLT.GetSize()):
         DYLT.SetBinError(j,(DYLT.GetBinError(j)*DYLT.GetBinError(j)+0.20*0.20*DYLT.GetBinContent(j)*DYLT.GetBinContent(j))**0.5)
      Data.Add(DYLT,-1.0)
      print DYLT.Integral(),Data.Integral()
      Data.SetName(dir_qcd[i])
      for k in range(1,Data.GetSize()-1):
        if Data.GetBinContent(k)<0:
   	  Data.SetBinContent(k,0)
      fileDataSub.cd()
      Data.Write()


