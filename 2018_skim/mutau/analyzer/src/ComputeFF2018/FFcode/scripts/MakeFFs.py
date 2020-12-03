#!/usr/bin/env python

import argparse
import os
import ComputeFF2018.FFcode.RawFF as RawFF
import ComputeFF2018.FFcode.FFmvisclosure as FFmvisclosure
import ComputeFF2018.FFcode.Closure as Closure
import ComputeFF2018.FFcode.FFOSSScorrection as FFOSSScorrection
import ComputeFF2018.FFcode.ControlPlots as ControlPlots
import ComputeFF2018.FFcode.Draw_raw as Draw_raw
import ComputeFF2018.FFcode.Draw_raw_tt as Draw_raw_tt
import ComputeFF2018.FFcode.Draw_closure_tt as Draw_closure_tt
import ComputeFF2018.FFcode.Draw_closure_et as Draw_closure_et
import ComputeFF2018.FFcode.Draw_control_tt as Draw_control_tt

parser = argparse.ArgumentParser(description="Master fake factor making script")
parser.add_argument('--channel','-c',nargs="?",choices=['mt','et','tt'],help="Which channel?",required=True)
parser.add_argument('--year','-y',nargs="?",choices=["2016","2017","2018"],help="Which year?",required=True)

args = parser.parse_args()

if not os.path.isdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/plots_"+args.channel+"_"+args.year):
    os.mkdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/plots_"+args.channel+"_"+args.year)

if not os.path.isdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/files_corrOSSSFF_"+args.channel):
    os.mkdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/files_corrOSSSFF_"+args.channel)
if not os.path.isdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/files_rawFF_"+args.channel):
    os.mkdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/files_rawFF_"+args.channel)
if not os.path.isdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/files_corr1FF_"+args.channel):
    os.mkdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/files_corr1FF_"+args.channel)
if args.channel=="tt":
  if not os.path.isdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/files_control_"+args.channel):
    os.mkdir(os.environ['CMSSW_BASE']+"/src/ComputeFF2018/files_control_"+args.channel)

os.system("rm "+os.environ["CMSSW_BASE"]+"/src/ComputeFF2018/files_corrOSSSFF_"+args.channel+"/*.root")
os.system("rm "+os.environ["CMSSW_BASE"]+"/src/ComputeFF2018/files_rawFF_"+args.channel+"/*.root")
os.system("rm "+os.environ["CMSSW_BASE"]+"/src/ComputeFF2018/files_corr1FF_"+args.channel+"/*.root")
os.system("rm "+os.environ["CMSSW_BASE"]+"/src/ComputeFF2018/FF_corrections_1.root")
os.system("rm "+os.environ["CMSSW_BASE"]+"/src/ComputeFF2018/FF_QCDcorrectionOSSS.root")
os.system("rm "+os.environ["CMSSW_BASE"]+"/src/ComputeFF2018/uncorrected_fakefactors_"+args.channel+".root")

if args.channel=="tt":
   os.system("rm "+os.environ["CMSSW_BASE"]+"/src/ComputeFF2018/files_control_"+args.channel+"/*.root")

RawFF.RawFF(args)
FFmvisclosure.FFmvisclosure(args)
Closure.Closure(args)
FFOSSScorrection.FFOSSScorrection(args)

if args.channel=="tt":
  ControlPlots.ControlPlots(args)
  Draw_raw_tt.Draw_raw_tt("raw",args.year)
  Draw_raw_tt.Draw_raw_tt("osss",args.year)
  Draw_closure_tt.Draw_closure_tt(args.year,"before")
  Draw_closure_tt.Draw_closure_tt(args.year,"after")
  Draw_control_tt.Draw_control_tt(args.year)

else:
  Draw_raw.Draw_raw("raw",args.year,args.channel)
  Draw_raw.Draw_raw("mvisclosure",args.year,args.channel)
  Draw_raw.Draw_raw("osss",args.year,args.channel)
  Draw_closure_et.Draw_closure_et(args.year,"before",args.channel)
  Draw_closure_et.Draw_closure_et(args.year,"after",args.channel)

os.system("mv *.pdf plots_"+args.channel+"_"+args.year)

finalPath = os.environ['CMSSW_BASE']+"/src/ComputeFF2018/ff_files_"+args.channel+"_"+args.year
if not os.path.isdir(finalPath):
    os.mkdir(finalPath)
os.system("cp uncorrected_fakefactors_"+args.channel+".root "+finalPath)
os.system("cp FF_corrections_1.root "+finalPath)
if args.channel=="et" or args.channel=="mt":
  os.system("cp FF_QCDcorrectionOSSS.root "+finalPath)
if args.channel=="tt":
  os.system("cp FF_QCDcorrectionOSSS_tt.root "+finalPath)
os.system("cp raw_FF_"+args.channel+".root "+finalPath)
os.system("cp mvisclosure_"+args.channel+".root "+finalPath)
os.system("cp OSSScorr_"+args.channel+".root "+finalPath)
