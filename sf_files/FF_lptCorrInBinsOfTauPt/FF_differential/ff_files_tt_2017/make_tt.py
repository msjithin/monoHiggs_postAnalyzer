import ROOT
import re
from array import array
import argparse

def make_tt():
   
   inFile = ROOT.TFile("raw_FF_tt.root","r")
   fileDataSub = ROOT.TFile("raw_FF_tautau_unfitted.root","RECREATE")
   dir_qcd=["tt_0jet_qcd_anti","tt_0jet_qcd_iso","tt_0SSloose_qcd_anti","tt_0SSloose_qcd_iso","tt_1jet_qcd_anti","tt_1jet_qcd_iso","tt_1SSloose_qcd_anti","tt_1SSloose_qcd_iso","tt_2jet_qcd_anti","tt_2jet_qcd_iso","tt_2SSloose_qcd_anti","tt_2SSloose_qcd_iso"]

   hisNames = ['raw_FF_tautau_0jet', 'raw_FF_tautau_1jet', 'raw_FF_tautau_2jet', 'raw_FF_tautau_3jet' ]
   ncat_qcd=12
      
   for i in range (0,11, 2):
      #print i, i+1
      print dir_qcd[i] , dir_qcd[i+1]
      Data_anti=inFile.Get(dir_qcd[i]).Get("data_obs")
      Data_anti.Add(inFile.Get(dir_qcd[i]).Get("VVLT"),-1.0)
      Data_anti.Add(inFile.Get(dir_qcd[i]).Get("VVJ"),-1.0)
      Data_anti.Add(inFile.Get(dir_qcd[i]).Get("TTJ"),-1.0)
      Data_anti.Add(inFile.Get(dir_qcd[i]).Get("TTLT"),-1.0)
      Data_anti.Add(inFile.Get(dir_qcd[i]).Get("DYJ"),-1.0)
      #Data_anti.Add(inFile.Get(dir_qcd[i]).Get("DYLT"), -1.0)
      DYLT=inFile.Get(dir_qcd[i]).Get("DYLT").Clone()
      for j in range(1,DYLT.GetSize()):
         DYLT.SetBinError(j,(DYLT.GetBinError(j)*DYLT.GetBinError(j)+0.20*0.20*DYLT.GetBinContent(j)*DYLT.GetBinContent(j))**0.5)
      Data_anti.Add(DYLT,-1.0)

      
      Data_iso=inFile.Get(dir_qcd[i+1]).Get("data_obs")
      Data_iso.Add(inFile.Get(dir_qcd[i+1]).Get("VVLT"),-1.0)
      Data_iso.Add(inFile.Get(dir_qcd[i+1]).Get("VVJ"),-1.0)
      Data_iso.Add(inFile.Get(dir_qcd[i+1]).Get("TTJ"),-1.0)
      Data_iso.Add(inFile.Get(dir_qcd[i+1]).Get("TTLT"),-1.0)
      Data_iso.Add(inFile.Get(dir_qcd[i+1]).Get("DYJ"),-1.0)
      Data_iso.Add(inFile.Get(dir_qcd[i+1]).Get("DYLT"), -1.0)
      DYLT=inFile.Get(dir_qcd[i+1]).Get("DYLT").Clone()
      for j in range(1,DYLT.GetSize()):
         DYLT.SetBinError(j,(DYLT.GetBinError(j)*DYLT.GetBinError(j)+0.20*0.20*DYLT.GetBinContent(j)*DYLT.GetBinContent(j))**0.5)
      Data_iso.Add(DYLT,-1.0)
      
      
      if '0jet' in dir_qcd[i]:
         histName = 'raw_FF_tautau_'+'0jet'
      elif '1jet' in dir_qcd[i]:
         histName = 'raw_FF_tautau_'+'1jet'
      elif '2jet' in dir_qcd[i]:
         histName = 'raw_FF_tautau_'+'2jet'
      else:
         histName = dir_qcd[i]

      Data_iso.Divide(Data_anti)
      Data_iso.SetName(histName)
      Data_iso.SetTitle(histName)
      fileDataSub.cd()
      Data_iso.Write()

make_tt()
  # KEY: TH1F     raw_FF_tautau_0jet;1    raw_FF_tautau_0jet
  # KEY: TH1F     raw_FF_tautau_1jet;1    raw_FF_tautau_1jet
  # KEY: TH1F     raw_FF_tautau_2jet;1    raw_FF_tautau_2jet
  # KEY: TH1F     raw_FF_tautau_3jet;1    raw_FF_tautau_3jet
  # KEY: TH1F     raw_FF_tautau_4jet;1    raw_FF_tautau_4jet
#  TFile*         raw_FF_tt.root
#   KEY: TDirectoryFile   tt_0jet_qcd_iso;1       tt_0jet_qcd_iso
#   KEY: TDirectoryFile   tt_0jet_qcd_anti;1      tt_0jet_qcd_anti
#   KEY: TDirectoryFile   tt_0jet_qcd_dm0_iso;1   tt_0jet_qcd_dm0_iso
#   KEY: TDirectoryFile   tt_0jet_qcd_dm0_anti;1  tt_0jet_qcd_dm0_anti
#   KEY: TDirectoryFile   tt_0jet_qcd_dm1_iso;1   tt_0jet_qcd_dm1_iso
#   KEY: TDirectoryFile   tt_0jet_qcd_dm1_anti;1  tt_0jet_qcd_dm1_anti
#   KEY: TDirectoryFile   tt_0jet_qcd_dm10_iso;1  tt_0jet_qcd_dm10_iso
#   KEY: TDirectoryFile   tt_0jet_qcd_dm10_anti;1 tt_0jet_qcd_dm10_anti
#   KEY: TDirectoryFile   tt_0jet_qcd_dm11_iso;1  tt_0jet_qcd_dm11_iso
#   KEY: TDirectoryFile   tt_0jet_qcd_dm11_anti;1 tt_0jet_qcd_dm11_anti
#   KEY: TDirectoryFile   tt_1jet_qcd_iso;1       tt_1jet_qcd_iso
#   KEY: TDirectoryFile   tt_1jet_qcd_anti;1      tt_1jet_qcd_anti
#   KEY: TDirectoryFile   tt_1jet_qcd_dm0_iso;1   tt_1jet_qcd_dm0_iso
#   KEY: TDirectoryFile   tt_1jet_qcd_dm0_anti;1  tt_1jet_qcd_dm0_anti
#   KEY: TDirectoryFile   tt_1jet_qcd_dm1_iso;1   tt_1jet_qcd_dm1_iso
#   KEY: TDirectoryFile   tt_1jet_qcd_dm1_anti;1  tt_1jet_qcd_dm1_anti
#   KEY: TDirectoryFile   tt_1jet_qcd_dm10_iso;1  tt_1jet_qcd_dm10_iso
#   KEY: TDirectoryFile   tt_1jet_qcd_dm10_anti;1 tt_1jet_qcd_dm10_anti
#   KEY: TDirectoryFile   tt_1jet_qcd_dm11_iso;1  tt_1jet_qcd_dm11_iso
#   KEY: TDirectoryFile   tt_1jet_qcd_dm11_anti;1 tt_1jet_qcd_dm11_anti
#   KEY: TDirectoryFile   tt_2jet_qcd_iso;1       tt_2jet_qcd_iso
#   KEY: TDirectoryFile   tt_2jet_qcd_anti;1      tt_2jet_qcd_anti
#   KEY: TDirectoryFile   tt_3jet_qcd_iso;1       tt_3jet_qcd_iso
#   KEY: TDirectoryFile   tt_3jet_qcd_anti;1      tt_3jet_qcd_anti
#   KEY: TDirectoryFile   tt_4jet_qcd_iso;1       tt_4jet_qcd_iso
#   KEY: TDirectoryFile   tt_4jet_qcd_anti;1      tt_4jet_qcd_anti
#   KEY: TDirectoryFile   tt_2jet_qcd_dm0_iso;1   tt_2jet_qcd_dm0_iso
#   KEY: TDirectoryFile   tt_2jet_qcd_dm0_anti;1  tt_2jet_qcd_dm0_anti
#   KEY: TDirectoryFile   tt_2jet_qcd_dm1_iso;1   tt_2jet_qcd_dm1_iso
#   KEY: TDirectoryFile   tt_2jet_qcd_dm1_anti;1  tt_2jet_qcd_dm1_anti
#   KEY: TDirectoryFile   tt_2jet_qcd_dm10_iso;1  tt_2jet_qcd_dm10_iso
#   KEY: TDirectoryFile   tt_2jet_qcd_dm10_anti;1 tt_2jet_qcd_dm10_anti
#   KEY: TDirectoryFile   tt_2jet_qcd_dm11_iso;1  tt_2jet_qcd_dm11_iso
#   KEY: TDirectoryFile   tt_2jet_qcd_dm11_anti;1 tt_2jet_qcd_dm11_anti
#   KEY: TDirectoryFile   tt_0SSloose_qcd_iso;1   tt_0SSloose_qcd_iso
#   KEY: TDirectoryFile   tt_0SSloose_qcd_anti;1  tt_0SSloose_qcd_anti
#   KEY: TDirectoryFile   tt_1SSloose_qcd_iso;1   tt_1SSloose_qcd_iso
#   KEY: TDirectoryFile   tt_1SSloose_qcd_anti;1  tt_1SSloose_qcd_anti
#   KEY: TDirectoryFile   tt_2SSloose_qcd_iso;1   tt_2SSloose_qcd_iso
#   KEY: TDirectoryFile   tt_2SSloose_qcd_anti;1  tt_2SSloose_qcd_anti
# root [2] tt_0jet_qcd_iso->cd()
# (bool) true
# root [3] .ls
# TDirectoryFile*         tt_0jet_qcd_iso tt_0jet_qcd_iso
#  KEY: TH1F      data_obs;1      h0LT_qcd_iso
#  KEY: TH1F      DYLT;1  h0LT_qcd_iso
#  KEY: TH1F      DYJ;1   h0J_qcd_iso
#  KEY: TH1F      TTLT;1  h0LT_qcd_iso
#  KEY: TH1F      TTJ;1   h0J_qcd_iso
#  KEY: TH1F      VVLT;1  h0LT_qcd_iso
#  KEY: TH1F      VVJ;1   h0J_qcd_iso
#  KEY: TH1F      STLT;1  h0LT_qcd_iso
#  KEY: TH1F      STJ;1   h0J_qcd_iso
