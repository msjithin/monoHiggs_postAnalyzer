#!/usr/bin/env python
import ROOT
import re
from array import array
import argparse


def Subtract_prompt(directory,channel):
  fileData=ROOT.TFile(directory+"/Data.root","r") #Data histogram
  fileVV=ROOT.TFile(directory+"/VV.root","r") #MC histogram with real leptons
  fileDY=ROOT.TFile(directory+"/DY.root","r") #MC histogram with real leptons
  fileTT=ROOT.TFile(directory+"/TT.root","r") #MC histogram with real leptons
  fileW=ROOT.TFile(directory+"/W.root","r") #MC histogram with real leptons
  fileDataSub=ROOT.TFile(directory+"/DataSub.root","recreate") #Data without real leptons

  dir_qcd=[channel+"_0jet_qcd_anti",channel+"_0jet_qcd_iso",channel+"_0SSloose_qcd_anti",channel+"_0SSloose_qcd_iso",channel+"_1jet_qcd_anti",channel+"_1jet_qcd_iso",channel+"_1SSloose_qcd_anti",channel+"_1SSloose_qcd_iso",channel+"_2jet_qcd_anti",channel+"_2jet_qcd_iso",channel+"_2SSloose_qcd_anti",channel+"_2SSloose_qcd_iso"]
  ncat_qcd=12

  if "corr1" in directory:
    print("corr1 subtraction")
    dir_qcd=[channel+"_0jet_qcd_anti",channel+"_0jet_qcd_iso",channel+"_0SSloose_qcd_anti",channel+"_0SSloose_qcd_iso",channel+"_0jet_qcd_mvis_anti",channel+"_0jet_qcd_mvis_iso",channel+"_1jet_qcd_mvis_anti",channel+"_1jet_qcd_mvis_iso",channel+"_2jet_qcd_mvis_anti",channel+"_2jet_qcd_mvis_iso"]
    ncat_qcd=10

  if "OSSSF" in directory:    
    print("OSSSF subtraction")
    ncat_qcd=2

  for i in range (0,ncat_qcd):    
    Data=fileData.Get(dir_qcd[i]).Get("data_obs")
    Data.Add(fileVV.Get(dir_qcd[i]).Get("VVLT"),-1.0)
    Data.Add(fileVV.Get(dir_qcd[i]).Get("VVJ"),-1.0)
    W=fileW.Get(dir_qcd[i]).Get("W").Clone()
    for j in range(1,W.GetSize()):
      W.SetBinError(j,(W.GetBinError(j)*W.GetBinError(j)+0.20*0.20*W.GetBinContent(j)*W.GetBinContent(j))**0.5)
    Data.Add(W,-1.0)
    #Data.Add(fileW.Get(dir_qcd[i]).Get("W"),-1.0)
    Data.Add(fileTT.Get(dir_qcd[i]).Get("TTJ"),-1.0)
    Data.Add(fileTT.Get(dir_qcd[i]).Get("TTLT"),-1.0)
    Data.Add(fileDY.Get(dir_qcd[i]).Get("DYJ"),-1.0)
    DYLT=fileDY.Get(dir_qcd[i]).Get("DYLT").Clone()
    for j in range(1,DYLT.GetSize()):
      DYLT.SetBinError(j,(DYLT.GetBinError(j)*DYLT.GetBinError(j)+0.10*0.10*DYLT.GetBinContent(j)*DYLT.GetBinContent(j))**0.5)
    Data.Add(DYLT,-1.0)
    #Data.Add(fileDY.Get(dir_qcd[i]).Get("DYLT"),-1.0)
    Data.SetName(dir_qcd[i])
    for k in range(1,Data.GetSize()-1):
      if Data.GetBinContent(k)<0:
        Data.SetBinContent(k,0)
      fileDataSub.cd()
    Data.Write()


  dir_w=[channel+"_0jet_w_anti",channel+"_0jet_w_iso",channel+"_1jet_w_anti",channel+"_1jet_w_iso",channel+"_2jet_w_anti",channel+"_2jet_w_iso"]
  ncat_w=6
  if "corr1" in directory:
    dir_w=[channel+"_0jet_w_anti",channel+"_0jet_w_iso",channel+"_0jet_w_mvis_anti",channel+"_0jet_w_mvis_iso",channel+"_1jet_w_mvis_anti",channel+"_1jet_w_mvis_iso",channel+"_2jet_w_mvis_anti",channel+"_2jet_w_mvis_iso"]
    ncat_w=8
  if "OSSSF" in directory:
    ncat_w=0

  for i in range (0,ncat_w):
    Data=fileData.Get(dir_w[i]).Get("data_obs")
    Data.Add(fileVV.Get(dir_w[i]).Get("VVLT"),-1.0)
    Data.Add(fileVV.Get(dir_w[i]).Get("VVJ"),-1.0)
    Data.Add(fileTT.Get(dir_w[i]).Get("TTJ"),-1.0)
    Data.Add(fileTT.Get(dir_w[i]).Get("TTLT"),-1.0)
    Data.Add(fileDY.Get(dir_w[i]).Get("DYJ"),-1.0)
    #Data.Add(fileDY.Get(dir_w[i]).Get("DYLT"),-1.0)
    DYLT=fileDY.Get(dir_w[i]).Get("DYLT").Clone()
    for j in range(1,DYLT.GetSize()):
      DYLT.SetBinError(j,(DYLT.GetBinError(j)*DYLT.GetBinError(j)+0.10*0.10*DYLT.GetBinContent(j)*DYLT.GetBinContent(j))**0.5)
    Data.Add(DYLT,-1.0)
    Data.SetName(dir_w[i])
    for k in range(1,Data.GetSize()-1):
      if Data.GetBinContent(k)<0:
        Data.SetBinContent(k,0)
      fileDataSub.cd()
    Data.Write()

  dir_tt=[channel+"_0jet_tt_anti",channel+"_0jet_tt_iso"]
  ncat_tt=2
  if "OSSSF" in directory:
    ncat_tt=0
  for i in range (0,ncat_tt):
    Data=fileData.Get(dir_tt[i]).Get("data_obs")
    Data.Add(fileVV.Get(dir_tt[i]).Get("VVLT"),-1.0)
    Data.Add(fileVV.Get(dir_tt[i]).Get("VVJ"),-1.0)
    W=fileW.Get(dir_tt[i]).Get("W").Clone()
    for j in range(1,W.GetSize()):
      W.SetBinError(j,(W.GetBinError(j)*W.GetBinError(j)+0.20*0.20*W.GetBinContent(j)*W.GetBinContent(j))**0.5)
    Data.Add(W,-1.0)
    #Data.Add(fileW.Get(dir_tt[i]).Get("W"),-1.0)
    Data.Add(fileTT.Get(dir_tt[i]).Get("TTLT"),-1.0)
    Data.Add(fileDY.Get(dir_tt[i]).Get("DYJ"),-1.0)
    DYLT=fileDY.Get(dir_tt[i]).Get("DYLT").Clone()
    for j in range(1,DYLT.GetSize()):
      DYLT.SetBinError(j,(DYLT.GetBinError(j)*DYLT.GetBinError(j)+0.10*0.10*DYLT.GetBinContent(j)*DYLT.GetBinContent(j))**0.5)
    Data.Add(DYLT,-1.0)
    #Data.Add(fileDY.Get(dir_tt[i]).Get("DYLT"),-1.0)
    Data.SetName(dir_tt[i])
    for k in range(1,Data.GetSize()-1):
      if Data.GetBinContent(k)<0:
        Data.SetBinContent(k,0)
      fileDataSub.cd()
    Data.Write()
