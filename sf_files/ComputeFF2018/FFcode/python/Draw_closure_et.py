#!/usr/bin/env python
import ROOT
import re
from array import array
import argparse

def add_lumi(datayear):
    lowX=0.58
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.06)
    lumi.SetTextFont (   42 )
    if datayear=="2018":
       lumi.AddText("2018, 59.7 fb^{-1} (13 TeV)")
    elif datayear=="2017":
       lumi.AddText("2017, 41.5 fb^{-1} (13 TeV)")
    elif datayear=="2016":
       lumi.AddText("2016, 35.9 fb^{-1} (13 TeV)")
    return lumi

def add_CMS():
    lowX=0.21
    lowY=0.70
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.05)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS Preliminary")
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
        output = ROOT.TLegend(0.45, 0.64, 0.92, 0.86, "", "brNDC")
        output.SetNColumns(2)
        output.SetLineWidth(0)
        output.SetLineStyle(0)
        output.SetFillStyle(0)
        output.SetBorderSize(0)
        output.SetTextFont(62)
        return output

def Draw_closure_et(year,corrections,channel):

   ROOT.gStyle.SetFrameLineWidth(3)
   ROOT.gStyle.SetLineWidth(3)
   ROOT.gStyle.SetOptStat(0)
   
   c=ROOT.TCanvas("canvas","",0,0,600,600)
   c.cd()
   
   
   adapt=ROOT.gROOT.GetColor(12)
   new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
   trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)
   
   
   categories=[channel+"_0jet_qcd_cat0jetlow_iso",channel+"_0jet_qcd_cat0jetlow_anti",channel+"_0jet_qcd_cat0jethigh_iso",channel+"_0jet_qcd_cat0jethigh_anti",channel+"_0jet_qcd_catboosted1_iso",channel+"_0jet_qcd_catboosted1_anti",channel+"_0jet_qcd_catboosted2_iso",channel+"_0jet_qcd_catboosted2_anti",channel+"_0jet_qcd_catvbflow_iso",channel+"_0jet_qcd_catvbflow_anti",channel+"_0jet_qcd_catvbfhigh_iso",channel+"_0jet_qcd_catvbfhigh_anti",channel+"_0jet_w_cat0jetlow_iso",channel+"_0jet_w_cat0jetlow_anti",channel+"_0jet_w_cat0jethigh_iso",channel+"_0jet_w_cat0jethigh_anti",channel+"_0jet_w_catboosted1_iso",channel+"_0jet_w_catboosted1_anti",channel+"_0jet_w_catboosted2_iso",channel+"_0jet_w_catboosted2_anti",channel+"_0jet_w_catvbflow_iso",channel+"_0jet_w_catvbflow_anti",channel+"_0jet_w_catvbfhigh_iso",channel+"_0jet_w_catvbfhigh_anti",channel+"_0jet_qcd_iso",channel+"_0jet_qcd_anti",channel+"_0jet_w_iso",channel+"_0jet_w_anti",channel+"_0jet_qcd_taupt_iso",channel+"_0jet_qcd_taupt_anti",channel+"_0jet_w_taupt_iso",channel+"_0jet_w_taupt_anti",channel+"_0jet_qcd_mjj_iso",channel+"_0jet_qcd_mjj_anti",channel+"_0jet_w_mjj_iso",channel+"_0jet_w_mjj_anti",channel+"_0jet_qcd_pth_iso",channel+"_0jet_qcd_pth_anti",channel+"_0jet_w_pth_iso",channel+"_0jet_w_pth_anti"]
   ncat=40
   file=ROOT.TFile("mvisclosure_"+channel+".root","r")
   if corrections=="after":
     file=ROOT.TFile("mvisclosure_"+channel+"_afterCorrections.root","r")
   
   for i in range (0,ncat/2):
      Data1=file.Get(categories[2*i]).Get("data_obs")
      TT1=file.Get(categories[2*i]).Get("TTLT")
      VV1=file.Get(categories[2*i]).Get("VVLT")
      DY1=file.Get(categories[2*i]).Get("DYLT")
      TTJ1=file.Get(categories[2*i]).Get("TTJ")
      W1=file.Get(categories[2*i]).Get("W")
      VVJ1=file.Get(categories[2*i]).Get("VVJ")
      DYJ1=file.Get(categories[2*i]).Get("DYJ")
   
      DataA1=file.Get(categories[2*i+1]).Get("data_obs")
      TTA1=file.Get(categories[2*i+1]).Get("TTLT")
      VVA1=file.Get(categories[2*i+1]).Get("VVLT")
      DYA1=file.Get(categories[2*i+1]).Get("DYLT")
      TTJA1=file.Get(categories[2*i+1]).Get("TTJ")
      WA1=file.Get(categories[2*i+1]).Get("W")
      VVJA1=file.Get(categories[2*i+1]).Get("VVJ")
      DYJA1=file.Get(categories[2*i+1]).Get("DYJ")

      print Data1.Integral(),DataA1.Integral()

      DataA=DataA1.Clone()
      TTA=TTA1.Clone()
      VVA=VVA1.Clone()
      DYA=DYA1.Clone()
      TTJA=TTJA1.Clone()
      WA=WA1.Clone()
      VVJA=VVJA1.Clone()
      DYJA=DYJA1.Clone()

      Data=Data1.Clone()
      TT=TT1.Clone()
      VV=VV1.Clone()
      DY=DY1.Clone()
      TTJ=TTJ1.Clone()
      W=W1.Clone()
      VVJ=VVJ1.Clone()
      DYJ=DYJ1.Clone()

      if "cat" in categories[2*i]:
            DataA3=DataA.Clone()
            TTA3=TTA.Clone()
            VVA3=VVA.Clone()
            DYA3=DYA.Clone()
            DYJA3=DYJA.Clone()
            TTJA3=TTJA.Clone()
            WA3=WA.Clone()
            VVJA3=VVJA.Clone()
            Data3=Data.Clone()
            TT3=TT.Clone()
            VV3=VV.Clone()
            DY3=DY.Clone()
            DYJ3=DYJ.Clone()
            TTJ3=TTJ.Clone()
            W3=W.Clone()
            VVJ3=VVJ.Clone()
            nn=DataA3.GetSize()
            nx=DataA3.GetXaxis().GetNbins()
            ny=DataA3.GetYaxis().GetNbins()
            Data=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            TT=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            VV=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            DY=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            TTJ=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            W=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            VVJ=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            DYJ=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            DataA=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            TTA=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            VVA=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            DYA=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            TTJA=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            WA=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            VVJA=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            DYJA=ROOT.TH1F("h1d","h1d",nx*ny,0,nx*ny)
            dir1=nx
            dir2=ny
            l=0
            for j in range(1,dir1+1):
                for k in range(1,dir2+1):
                   l=l+1
                   n = DataA3.GetBin(j,k);
                   Data.SetBinContent(l,Data3.GetBinContent(n))
                   Data.SetBinError(l,Data3.GetBinError(n))
                   DY.SetBinContent(l,DY3.GetBinContent(n))
                   DY.SetBinError(l,DY3.GetBinError(n))
                   TT.SetBinContent(l,TT3.GetBinContent(n))
                   TT.SetBinError(l,TT3.GetBinError(n))
                   VV.SetBinContent(l,VV3.GetBinContent(n))
                   VV.SetBinError(l,VV3.GetBinError(n))
                   DYJ.SetBinContent(l,DYJ3.GetBinContent(n))
                   DYJ.SetBinError(l,DYJ3.GetBinError(n))
                   TTJ.SetBinContent(l,TTJ3.GetBinContent(n))
                   TTJ.SetBinError(l,TTJ3.GetBinError(n))
                   W.SetBinContent(l,W3.GetBinContent(n))
                   W.SetBinError(l,W3.GetBinError(n))
                   VVJ.SetBinContent(l,VVJ3.GetBinContent(n))
                   VVJ.SetBinError(l,VVJ3.GetBinError(n))

                   DataA.SetBinContent(l,DataA3.GetBinContent(n))
                   DataA.SetBinError(l,DataA3.GetBinError(n))
                   DYA.SetBinContent(l,DYA3.GetBinContent(n))
                   DYA.SetBinError(l,DYA3.GetBinError(n))
                   TTA.SetBinContent(l,TTA3.GetBinContent(n))
                   TTA.SetBinError(l,TTA3.GetBinError(n))
                   VVA.SetBinContent(l,VVA3.GetBinContent(n))
                   VVA.SetBinError(l,VVA3.GetBinError(n))
                   DYJA.SetBinContent(l,DYJA3.GetBinContent(n))
                   DYJA.SetBinError(l,DYJA3.GetBinError(n))
                   TTJA.SetBinContent(l,TTJA3.GetBinContent(n))
                   TTJA.SetBinError(l,TTJA3.GetBinError(n))
                   WA.SetBinContent(l,WA3.GetBinContent(n))
                   WA.SetBinError(l,WA3.GetBinError(n))
                   VVJA.SetBinContent(l,VVJA3.GetBinContent(n))
                   VVJA.SetBinError(l,VVJA3.GetBinError(n))

      if "_qcd_" in categories[2*i]:
	 print WA.Integral(),DYJA.Integral()
	 DataA.Add(VVA,-1)
         DataA.Add(VVJA,-1)
         for j in range(1,WA.GetSize()):
           WA.SetBinError(j,(WA.GetBinError(j)*WA.GetBinError(j)+0.20*0.20*WA.GetBinContent(j)*WA.GetBinContent(j))**0.5)
         DataA.Add(WA,-1.0)
         DataA.Add(TTA,-1)
         DataA.Add(TTJA,-1)
         DataA.Add(DYJA,-1)
         for j in range(1,DYA.GetSize()):
           DYA.SetBinError(j,(DYA.GetBinError(j)*DYA.GetBinError(j)+0.10*0.10*DYA.GetBinContent(j)*DYA.GetBinContent(j))**0.5)
         DataA.Add(DYA,-1.0)

         Data.Add(VV,-1)
         Data.Add(VVJ,-1)
         for j in range(1,W.GetSize()):
           W.SetBinError(j,(W.GetBinError(j)*W.GetBinError(j)+0.20*0.20*W.GetBinContent(j)*W.GetBinContent(j))**0.5)
         Data.Add(W,-1.0)
         Data.Add(TT,-1)
         Data.Add(TTJ,-1)
         Data.Add(DYJ,-1)
         for j in range(1,DY.GetSize()):
           DY.SetBinError(j,(DY.GetBinError(j)*DY.GetBinError(j)+0.10*0.10*DY.GetBinContent(j)*DY.GetBinContent(j))**0.5)
         Data.Add(DY,-1.0)

      if "_w_" in categories[2*i]:
         Data.Add(VV,-1)
         Data.Add(VVJ,-1)
         Data.Add(TT,-1)
         Data.Add(TTJ,-1)
         Data.Add(DYJ,-1)
         for j in range(1,DY.GetSize()):
           DY.SetBinError(j,(DY.GetBinError(j)*DY.GetBinError(j)+0.10*0.10*DY.GetBinContent(j)*DY.GetBinContent(j))**0.5)
         Data.Add(DY,-1.0)

         DataA.Add(VVA,-1)
         DataA.Add(VVJA,-1)
         DataA.Add(TTA,-1)
         DataA.Add(TTJA,-1)
         DataA.Add(DYJA,-1)
         for j in range(1,DYA.GetSize()):
           DYA.SetBinError(j,(DYA.GetBinError(j)*DYA.GetBinError(j)+0.10*0.10*DYA.GetBinContent(j)*DYA.GetBinContent(j))**0.5)
         DataA.Add(DYA,-1.0)
   
      Data.GetXaxis().SetTitle("")
      Data.GetXaxis().SetTitleSize(0)
      Data.GetXaxis().SetNdivisions(505)
      Data.GetYaxis().SetLabelFont(42)
      Data.GetYaxis().SetLabelOffset(0.01)
      Data.GetYaxis().SetLabelSize(0.06)
      Data.GetYaxis().SetTitleSize(0.075)
      Data.GetYaxis().SetTitleOffset(1.18)
      Data.SetTitle("")
      Data.GetYaxis().SetTitle("Events/bin")
   
      DataA.SetFillColor(ROOT.TColor.GetColor("#12cadd"))
      DataA.SetLineColor(1)
   
      Data.SetMarkerStyle(20)
      Data.SetMarkerSize(1)
      Data.SetLineColor(1)
      Data.SetLineWidth(2)
   
      errorBand = DataA.Clone()
      errorBand.SetMarkerSize(0)
      errorBand.SetFillColor(new_idx)
      errorBand.SetFillStyle(3001)
      errorBand.SetLineWidth(1)
   
      pad1 = ROOT.TPad("pad1","pad1",0,0.35,1,1)
      pad1.Draw()
      pad1.cd()
      pad1.SetFillColor(0)
      pad1.SetBorderMode(0)
      pad1.SetBorderSize(10)
      pad1.SetTickx(1)
      pad1.SetTicky(1)
      pad1.SetLeftMargin(0.18)
      pad1.SetRightMargin(0.05)
      pad1.SetTopMargin(0.122)
      pad1.SetBottomMargin(0.026)
      pad1.SetFrameFillStyle(0)
      pad1.SetFrameLineStyle(0)
      pad1.SetFrameLineWidth(3)
      pad1.SetFrameBorderMode(0)
      pad1.SetFrameBorderSize(10)
   
      Data.GetXaxis().SetLabelSize(0)
      Data.SetMaximum(max(Data.GetMaximum()*1.5,errorBand.GetMaximum()*1.5))
      Data.SetMinimum(0.0)
      Data.Draw("e0")
      DataA.Draw("histsame")
      errorBand.Draw("e2same")
      Data.Draw("e0same")

      categ2  = ROOT.TPaveText(0.65, 0.6+0.013, 0.93, 0.6+0.155, "NDC")
      categ2.SetBorderSize(   0 )
      categ2.SetFillStyle(    0 )
      categ2.SetTextAlign(   12 )
      categ2.SetTextSize ( 0.05 )
      categ2.SetTextColor(    1 )
      categ2.SetTextFont (   42 )
      categ2.AddText("KS prob. = %.5f" %(Data.KolmogorovTest(DataA)))
      categ2.Draw("same")

      categ3  = ROOT.TPaveText(0.65, 0.5+0.013, 0.93, 0.5+0.12, "NDC")
      categ3.SetBorderSize(   0 )
      categ3.SetFillStyle(    0 )
      categ3.SetTextAlign(   12 )
      categ3.SetTextSize ( 0.05 )
      categ3.SetTextColor(    1 )
      categ3.SetTextFont (   42 )
      categ3.AddText("Observed: %.1f" %(Data.Integral()))
      categ3.AddText("Predicted: %.1f" %(DataA.Integral()))
      print Data.Integral(),DataA.Integral()
      categ3.Draw("same")

   
      legende=make_legend()
      legende.AddEntry(Data,"Observed","elp")
      legende.AddEntry(DataA,"Predicted","f")
      legende.AddEntry(errorBand,"Stat. unc.","f")
      legende.Draw()
   
      l1=add_lumi(year)
      l1.Draw("same")
      l2=add_CMS()
      l2.Draw("same")
      l3=add_Preliminary()
    
      pad1.RedrawAxis()
   
      c.cd()
      pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.35);
      pad2.SetTopMargin(0.05);
      pad2.SetBottomMargin(0.35);
      pad2.SetLeftMargin(0.18);
      pad2.SetRightMargin(0.05);
      pad2.SetTickx(1)
      pad2.SetTicky(1)
      pad2.SetFrameLineWidth(3)
      pad2.SetGridx()
      pad2.SetGridy()
      pad2.Draw()
      pad2.cd()
      h1=Data.Clone()
      h1.SetMaximum(1.5)#FIXME(1.5)
      h1.SetMinimum(0.5)#FIXME(0.5)
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
      h1.GetXaxis().SetTitle("Bin number")
      if i==12 or i==13: 
        h1.GetXaxis().SetTitle("m_{vis} (GeV)")
      if i==14 or i==15:
        h1.GetXaxis().SetTitle("#tau_{h} p_{T} (GeV)")
      if i==16 or i==17:
        h1.GetXaxis().SetTitle("m_{jj} (GeV)")
      if i==18 or i==19:
        h1.GetXaxis().SetTitle("Higgs p_{T} (GeV)")

      h1.GetXaxis().SetLabelSize(0.08)
      h1.GetYaxis().SetLabelSize(0.08)
      h1.GetYaxis().SetTitle("Obs./Exp.")
      h1.GetXaxis().SetNdivisions(505)
      h1.GetYaxis().SetNdivisions(5)
   
      h1.GetXaxis().SetTitleSize(0.15)
      h1.GetYaxis().SetTitleSize(0.15)
      h1.GetYaxis().SetTitleOffset(0.56)
      h1.GetXaxis().SetTitleOffset(1.04)
      h1.GetXaxis().SetLabelSize(0.11)
      h1.GetYaxis().SetLabelSize(0.11)
      h1.GetXaxis().SetTitleFont(42)
      h1.GetYaxis().SetTitleFont(42)
   
      h1.Draw("e0p")
      h3.Draw("e2same")
   
      c.cd()
      pad1.Draw()
   
      ROOT.gPad.RedrawAxis()
   
      c.Modified()
      if corrections=="after":
         c.SaveAs("closure_compa_"+categories[2*i]+"_afterCorrections.pdf")
      else:
         c.SaveAs("closure_compa_"+categories[2*i]+".pdf")
   
