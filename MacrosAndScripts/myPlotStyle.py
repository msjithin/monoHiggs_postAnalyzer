import ROOT


def add_lumi(year, channel_):
  lowX=0.40
  lowY=0.825
  lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.50, lowY+0.16, "NDC")
  lumi.SetBorderSize(   0 )
  lumi.SetFillStyle(    0 )
  lumi.SetTextAlign(   32 )#12
  lumi.SetTextColor(    1 )
  lumi.SetTextSize(0.05)
  lumi.SetTextFont (   42 )
  lumiProcessed="41.52"
  if year=="2018":
    lumiProcessed="59.7"
  if channel_=="combined":
    lumi.AddText("4 channels combined "+year+" , "+lumiProcessed+" fb^{-1} (13 TeV)")
  if channel_=="mutau":
    lumi.AddText("#mu#tau_{h} "+year+" , "+lumiProcessed+" fb^{-1} (13 TeV)")
  if channel_=="etau":
    lumi.AddText("e#tau_{h} "+year+" , "+lumiProcessed+" fb^{-1} (13 TeV)")
  if channel_=="tautau":
    lumi.AddText("#tau_{h}#tau_{h} "+year+" , "+lumiProcessed+" fb^{-1} (13 TeV)")
  if channel_=="emu":
    lumi.AddText("e#mu "+year+" , "+lumiProcessed+" fb^{-1} (13 TeV)")
  return lumi

def add_CMS():
    lowX=0.65 #0.21
    lowY=0.68
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.08)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   31 )#12
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
    return lumi

def add_Preliminary():
    lowX=0.65 # 0.21
    lowY=0.63
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextFont (   40 )
    lumi.SetTextSize(0.06)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   31 )#12
    lumi.SetTextColor(    1 )
    lumi.AddText("Preliminary")
    return lumi

def make_legend():
  output = ROOT.TLegend(0.85, 0.45, 1.0, 0.75, "", "brNDC")
  #output = ROOT.TLegend(0.2, 0.1, 0.47, 0.65, "", "brNDC")
  output.SetLineWidth(1)
  output.SetLineStyle(1)
  output.SetFillStyle(1001) #0
  output.SetFillColor(0)
  output.SetBorderSize(1)
  output.SetTextFont(42)
  return output


# declare colors
color_ztt="#ffcc66"
color_zll="#4496c8"
color_tt="#9999cc"
color_ggh="#12cadd"
color_vv="#990099"
color_wjets="#cc6666"


errorStyle=3002

c=ROOT.TCanvas("canvas","",0,0,1300,1200)
pad1 = ROOT.TPad("pad1","pad1",0,0.25,1,1)
pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetBorderSize(1)
pad1.SetTickx(1)
pad1.SetTicky(1)
pad1.SetGridx()
pad1.SetLeftMargin(0.15) #0.15
pad1.SetRightMargin(0.15) #0.1
pad1.SetTopMargin(0.122)
pad1.SetBottomMargin(0.025)
pad1.SetFrameFillStyle(0)
pad1.SetFrameLineStyle(0)
pad1.SetFrameLineWidth(1)
pad1.SetFrameBorderMode(0)
pad1.SetFrameBorderSize(1)

categ  = ROOT.TPaveText(0.21, 0.5+0.013, 0.43, 0.70+0.155, "NDC")
categ.SetBorderSize(   0 )
categ.SetFillStyle(    0 )
categ.SetTextAlign(   12 )
categ.SetTextSize ( 0.06 )
categ.SetTextColor(    1 )
categ.SetTextFont (   42 )

pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.25);
pad2.SetTopMargin(0.02);
pad2.SetBottomMargin(0.35);
pad2.SetLeftMargin(0.15);
pad2.SetRightMargin(0.15);
pad2.SetTickx(1)
pad2.SetTicky(1)
pad2.SetFrameLineWidth(1)
#pad2.SetGridx()
pad2.SetGridy()
#pad2.SetLogy()
