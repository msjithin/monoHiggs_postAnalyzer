#!/usr/bin/env python
import argparse
import os
import ComputeFF2018.FFcode.Subtract_prompt as Subtract_prompt
import ComputeFF2018.FFcode.Subtract_prompt_tt as Subtract_prompt_tt

def FFmvisclosure(args):
    if args.year == "2016":
        if args.channel == "mt":
            path = '/data/ccaillol/smhmt2016_svfitted_12oct/'
        if args.channel == "et":
            path = '/data/ccaillol/smhet2016_svfitted_12oct/'
        if args.channel == "tt":
            path = '/data/ccaillol/smhtt2016_12oct/'
    elif args.year == "2017":
        if args.channel == "mt":
            path = '/data/ccaillol/smhmt2017_svfitted_12oct/'
        if args.channel == "et":
            path = '/data/ccaillol/smhet2017_svfitted_12oct/'
        if args.channel == "tt":
            path = '/data/ccaillol/smhtt2017_12oct/'
    elif args.year == "2018":
        if args.channel == "mt":
            path = '/data/ccaillol/smhmt2018_svfitted_12oct/'
        if args.channel == "et":
            path = '/data/ccaillol/smhet2018_svfitted_12oct/'
        if args.channel == "tt":
            path = '/data/ccaillol/smhtt2018_12oct/'

    #figure out our executable and our output directory
    if args.channel == 'mt':
        executable = "Set1_correction_mt"
        outputPath = os.environ['CMSSW_BASE']+'/src/ComputeFF2018/files_corr1FF_mt/'
    elif args.channel == "et":
        executable = "Set1_correction_et"
        outputPath = os.environ['CMSSW_BASE']+'/src/ComputeFF2018/files_corr1FF_et/'
    elif args.channel == "tt":
        executable = "Set1_correction_tt"
        outputPath = os.environ['CMSSW_BASE']+'/src/ComputeFF2018/files_corr1FF_tt/'

    if args.year == '2018':
        commandParams = [
            [executable,path+'DataA.root',outputPath+'DataA.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataB.root',outputPath+'DataB.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataC.root',outputPath+'DataC.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataD.root',outputPath+'DataD.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DY.root',outputPath+'DYincl.root','DY','DY',args.year,'no'],
            [executable,path+'DY1.root',outputPath+'DY1.root','DY','DY',args.year,'no'],
            [executable,path+'DY2.root',outputPath+'DY2.root','DY','DY',args.year,'no'],
            [executable,path+'DY3.root',outputPath+'DY3.root','DY','DY',args.year,'no'],
            [executable,path+'DY4.root',outputPath+'DY4.root','DY','DY',args.year,'no'],
            [executable,path+'TTToHadronic.root',outputPath+'TTToHadronic.root','TTToHadronic','TT',args.year,'no'],
            [executable,path+'TTTo2L2Nu.root',outputPath+'TTTo2L2Nu.root',' TTTo2L2Nu','TT',args.year,'no'],
            [executable,path+'TTToSemiLeptonic.root',outputPath+'TTToSemiLeptonic.root','TTToSemiLeptonic','TT',args.year,'no'],
            [executable,path+'WW.root',outputPath+'WW.root','WW','VV',args.year,'no'],
            [executable,path+'WZ.root',outputPath+'WZ.root','WZ','VV',args.year,'no'],
            [executable,path+'ZZ.root',outputPath+'ZZ.root','ZZ','VV',args.year,'no'],
            [executable,path+'ST_t_antitop.root',outputPath+'ST_t_antitop.root','ST_t_antitop','ST',args.year,'no'],
            [executable,path+'ST_t_top.root',outputPath+'ST_t_top.root','ST_t_top','ST',args.year,'no'],
            [executable,path+'ST_tW_antitop.root',outputPath+'ST_tW_antitop.root','ST_tW_antitop','ST',args.year,'no'],
            [executable,path+'ST_tW_top.root',outputPath+'ST_tW_top.root','ST_tW_top','ST',args.year,'no'],
        ]
        if args.channel == "et" or args.channel=="mt":
            commandParams.append([executable,path+'Wall.root',outputPath+'W.root','W','W',args.year,'no'])
            commandParams.append([executable,path+'Wall.root',outputPath+'WMC.root','W','WMC',args.year,'no'])
            commandParams.append([executable,path+'Wall.root',outputPath+'WMC2.root','W','WMC2',args.year,'no'])
            commandParams.append([executable,path+'TTToHadronic.root',outputPath+'TTToHadronicMC.root','TTToHadronic','TTMC',args.year,'no'])
            commandParams.append([executable,path+'TTTo2L2Nu.root',outputPath+'TTTo2L2NuMC.root',' TTTo2L2Nu','TTMC',args.year,'no'])
            commandParams.append([executable,path+'TTToSemiLeptonic.root',outputPath+'TTToSemiLeptonicMC.root','TTToSemiLeptonic','TTMC',args.year,'no'])
        haddFiles = {
            "Data.root": [outputPath+"DataA.root",outputPath+"DataB.root",outputPath+"DataC.root",outputPath+"DataD.root"],
            "DY.root": [outputPath+"DYincl.root",outputPath+"DY1.root",outputPath+"DY2.root",outputPath+"DY3.root",outputPath+"DY4.root"],
            "TT.root": [outputPath+"TTToHadronic.root",outputPath+"TTToSemiLeptonic.root",outputPath+"TTTo2L2Nu.root"],
            "VV.root": [outputPath+"WW.root",outputPath+"WZ.root",outputPath+"ZZ.root",outputPath+"ST_t_antitop.root",outputPath+"ST_t_top.root",outputPath+"ST_tW_antitop.root",outputPath+"ST_tW_top.root"]
        }
	if args.channel == "et" or args.channel=="mt":
           haddFiles = {
               "Data.root": [outputPath+"DataA.root",outputPath+"DataB.root",outputPath+"DataC.root",outputPath+"DataD.root"],
               "DY.root": [outputPath+"DYincl.root",outputPath+"DY1.root",outputPath+"DY2.root",outputPath+"DY3.root",outputPath+"DY4.root"],
               "TT.root": [outputPath+"TTToHadronic.root",outputPath+"TTToSemiLeptonic.root",outputPath+"TTTo2L2Nu.root"],
               "VV.root": [outputPath+"WW.root",outputPath+"WZ.root",outputPath+"ZZ.root",outputPath+"ST_t_antitop.root",outputPath+"ST_t_top.root",outputPath+"ST_tW_antitop.root",outputPath+"ST_tW_top.root"],
	       "TTMC.root": [outputPath+"TTToHadronicMC.root",outputPath+"TTToSemiLeptonicMC.root",outputPath+"TTTo2L2NuMC.root"]
           }
    elif args.year=='2017':
        commandParams = [
            [executable,path+'DataB.root',outputPath+'DataB.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataC.root',outputPath+'DataC.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataD.root',outputPath+'DataD.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataE.root',outputPath+'DataE.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataF.root',outputPath+'DataF.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DY.root',outputPath+'DYincl.root','DY','DY',args.year,'no'],
            [executable,path+'DY1.root',outputPath+'DY1.root','DY','DY',args.year,'no'],
            [executable,path+'DY2.root',outputPath+'DY2.root','DY','DY',args.year,'no'],
            [executable,path+'DY3.root',outputPath+'DY3.root','DY','DY',args.year,'no'],
            [executable,path+'DY4.root',outputPath+'DY4.root','DY','DY',args.year,'no'],
            [executable,path+'TTToHadronic.root',outputPath+'TTToHadronic.root','TTToHadronic','TT',args.year,'no'],
            [executable,path+'TTTo2L2Nu.root',outputPath+'TTTo2L2Nu.root',' TTTo2L2Nu','TT',args.year,'no'],
            [executable,path+'TTToSemiLeptonic.root',outputPath+'TTToSemiLeptonic.root','TTToSemiLeptonic','TT',args.year,'no'],
            [executable,path+'WW.root',outputPath+'WW.root','WW','VV',args.year,'no'],
            [executable,path+'WZ.root',outputPath+'WZ.root','WZ','VV',args.year,'no'],
            [executable,path+'ZZ.root',outputPath+'ZZ.root','ZZ','VV',args.year,'no'],
            [executable,path+'ST_t_antitop.root',outputPath+'ST_t_antitop.root','ST_t_antitop','ST',args.year,'no'],
            [executable,path+'ST_t_top.root',outputPath+'ST_t_top.root','ST_t_top','ST',args.year,'no'],
            [executable,path+'ST_tW_antitop.root',outputPath+'ST_tW_antitop.root','ST_tW_antitop','ST',args.year,'no'],
            [executable,path+'ST_tW_top.root',outputPath+'ST_tW_top.root','ST_tW_top','ST',args.year,'no']
        ]
        if args.channel == "et" or args.channel=="mt":
            commandParams.append([executable,path+'Wall.root',outputPath+'W.root','W','W',args.year,'no'])
            commandParams.append([executable,path+'Wall.root',outputPath+'WMC.root','W','WMC',args.year,'no'])
            commandParams.append([executable,path+'Wall.root',outputPath+'WMC2.root','W','WMC2',args.year,'no'])
            commandParams.append([executable,path+'TTToHadronic.root',outputPath+'TTToHadronicMC.root','TTToHadronic','TTMC',args.year,'no'])
            commandParams.append([executable,path+'TTTo2L2Nu.root',outputPath+'TTTo2L2NuMC.root',' TTTo2L2Nu','TTMC',args.year,'no'])
            commandParams.append([executable,path+'TTToSemiLeptonic.root',outputPath+'TTToSemiLeptonicMC.root','TTToSemiLeptonic','TTMC',args.year,'no'])
        haddFiles = {
            "Data.root": [outputPath+"DataB.root",outputPath+"DataC.root",outputPath+"DataD.root",outputPath+"DataE.root",outputPath+"DataF.root"],
            "DY.root": [outputPath+"DYincl.root",outputPath+"DY1.root",outputPath+"DY2.root",outputPath+"DY3.root",outputPath+"DY4.root"],
            "TT.root": [outputPath+"TTToHadronic.root",outputPath+"TTToSemiLeptonic.root",outputPath+"TTTo2L2Nu.root"],
            "VV.root": [outputPath+"WW.root",outputPath+"WZ.root",outputPath+"ZZ.root",outputPath+"ST_t_antitop.root",outputPath+"ST_t_top.root",outputPath+"ST_tW_antitop.root",outputPath+"ST_tW_top.root"]
        }
        if args.channel == "et" or args.channel=="mt":
          haddFiles = {
              "Data.root": [outputPath+"DataB.root",outputPath+"DataC.root",outputPath+"DataD.root",outputPath+"DataE.root",outputPath+"DataF.root"],
              "DY.root": [outputPath+"DYincl.root",outputPath+"DY1.root",outputPath+"DY2.root",outputPath+"DY3.root",outputPath+"DY4.root"],
              "TT.root": [outputPath+"TTToHadronic.root",outputPath+"TTToSemiLeptonic.root",outputPath+"TTTo2L2Nu.root"],
              "VV.root": [outputPath+"WW.root",outputPath+"WZ.root",outputPath+"ZZ.root",outputPath+"ST_t_antitop.root",outputPath+"ST_t_top.root",outputPath+"ST_tW_antitop.root",outputPath+"ST_tW_top.root"],
	      "TTMC.root": [outputPath+"TTToHadronicMC.root",outputPath+"TTToSemiLeptonicMC.root",outputPath+"TTTo2L2NuMC.root"]
          }
    elif args.year=='2016':
        commandParams = [
            [executable,path+'DataB.root',outputPath+'DataB.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataC.root',outputPath+'DataC.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataD.root',outputPath+'DataD.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataE.root',outputPath+'DataE.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataF.root',outputPath+'DataF.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataG.root',outputPath+'DataG.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DataH.root',outputPath+'DataH.root','data_obs','data_obs',args.year,'no'],
            [executable,path+'DY.root',outputPath+'DYincl.root','DY','DY',args.year,'no'],
            [executable,path+'DY1.root',outputPath+'DY1.root','DY','DY',args.year,'no'],
            [executable,path+'DY2.root',outputPath+'DY2.root','DY','DY',args.year,'no'],
            [executable,path+'DY3.root',outputPath+'DY3.root','DY','DY',args.year,'no'],
            [executable,path+'DY4.root',outputPath+'DY4.root','DY','DY',args.year,'no'],
            [executable,path+'TT.root',outputPath+'TT.root','TT','TT',args.year,'no'],
            [executable,path+'WW.root',outputPath+'WW.root','WW','VV',args.year,'no'],
            [executable,path+'WZ.root',outputPath+'WZ.root','WZ','VV',args.year,'no'],
            [executable,path+'ZZ.root',outputPath+'ZZ.root','ZZ','VV',args.year,'no'],
            [executable,path+'ST_t_antitop.root',outputPath+'ST_t_antitop.root','ST_t_antitop','ST',args.year,'no'],
            [executable,path+'ST_t_top.root',outputPath+'ST_t_top.root','ST_t_top','ST',args.year,'no'],
            [executable,path+'ST_tW_antitop.root',outputPath+'ST_tW_antitop.root','ST_tW_antitop','ST',args.year,'no'],
            [executable,path+'ST_tW_top.root',outputPath+'ST_tW_top.root','ST_tW_top','ST',args.year,'no'],            
            [executable,path+'ggH125.root',outputPath+'ggH_htt125.root','ggH_htt125','ggH_htt125',args.year,'no'],
        ]
	if args.channel=="mt" or args.channel=="et":
            commandParams.append([executable,path+'Wall.root',outputPath+'W.root','W','W',args.year,'no'])
            commandParams.append([executable,path+'Wall.root',outputPath+'WMC.root','W','WMC',args.year,'no'])
            commandParams.append([executable,path+'TT.root',outputPath+'TTMC.root','TT','TTMC',args.year,'no'])
        haddFiles = {
            "Data.root": [outputPath+"DataB.root",outputPath+"DataC.root",outputPath+"DataD.root",outputPath+"DataE.root",outputPath+"DataF.root",outputPath+"DataG.root",outputPath+"DataH.root"],
            "DY.root": [outputPath+"DYincl.root",outputPath+"DY1.root",outputPath+"DY2.root",outputPath+"DY3.root",outputPath+"DY4.root"],
            "VV.root": [outputPath+"WW.root",outputPath+"WZ.root",outputPath+"ZZ.root",outputPath+"ST_t_antitop.root",outputPath+"ST_t_top.root",outputPath+"ST_tW_antitop.root",outputPath+"ST_tW_top.root"]
        }
    #Run all sets
    for command in commandParams:
        commandString = ''
        for element in command:
            commandString+=element+' '
        os.system(commandString)
    #hadd them
    for haddedFile in haddFiles:
        haddCommand = "hadd -f "+outputPath+haddedFile+' '
        for inputFile in haddFiles[haddedFile]:
            haddCommand+=inputFile+' '
        os.system(haddCommand)
    #do our subtractions, and our fitting.
    if args.channel=="et" or args.channel=="mt":
      Subtract_prompt.Subtract_prompt(outputPath,args.channel)
    elif args.channel=="tt":
      Subtract_prompt_tt.Subtract_prompt_tt(outputPath)
    os.system("root -l -b -q \'Fit_FFclosure_"+args.channel+".cc("+args.year+")\'")
    
    #this one has to specifically be done here because the file for weighting doesn't even exist until the fits are done.
    if (args.channel=="et" or args.channel=="mt"):
      os.system(executable+" "+path+"Wall.root"+" "+outputPath+"WMC2.root W WMC2 "+args.year+" no")
    
    #make the final mvisclosure file
    finalMvisClosureCommand = "hadd -f mvisclosure_"+args.channel+".root "
    finalMvisClosureCommand+=outputPath+"Data.root "
    finalMvisClosureCommand+=outputPath+"DY.root "
    if args.channel=="et" or args.channel=="mt":
      finalMvisClosureCommand+=outputPath+"W.root "
    finalMvisClosureCommand+=outputPath+"TT.root "
    finalMvisClosureCommand+=outputPath+"VV.root "
    os.system(finalMvisClosureCommand)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Master fake factor making script")
    parser.add_argument('--channel','-c',nargs="?",choices=['mt','et','tt'],help="Which channel?",required=True)
    parser.add_argument('--year','-y',nargs="?",choices=["2016","2017","2018"],help="Which year?",required=True)

    args=parser.parse_args()
    FFmvisclosure(args)
    
