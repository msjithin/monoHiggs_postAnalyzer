#!/usr/bin/python
import subprocess
from os import path, listdir, getcwd, kill, popen
import argparse
import signal
import time 

current_path = getcwd()


var_mapping = { 1: "elePt" ,     2:"eleEta",          3:"elePhi" ,
                4: "muPt" ,      5:"muEta",           6:"muPhi" ,
                7:"tauPt",       8:"tauEta",          9:"tauPhi",
                10:"tau1Pt",     11:"tau1Eta",        12:"tau1Phi",
                13:"tau2Pt",     14:"tau2Eta",        15:"tau2Phi",
                16:"tauIso",     17:"tauDecayMode",   18:"tauCharge", 
                19:"deltaR",     20:"higgsPt",        21:"nJet", 
                22:"visMass",    23:"mT_eleMet",      24:"trigger",  
                25:"genMatch",   26:"met",            27:"metPhi", 
                28:"deltaPhi",   29:"deltaEta",       30:"metLongXaxis",  
                31:"tot_TMass",  32:"tot_TMass_full", 33 : "metFull", 
                34 :"etau",      35:"mutau",          36:"tautau",
                41: "muPt_new",  42 : "tauPt_new", 
                43: "muMass_new",44: "tauMass_new", 
                45:"met_new",    46:"visMass_new",
                47:"higgsPt_new" , 48:"visMass_met_new", 49:"higgsPt_met_new",
                50: "muMetTMassFull",
                }

title_mapping= {
                "tot_TMass" : "tot tr mass",  
                "tot_TMass_full" : "tot tr mass", 
                "metFull" : "MET", 
}
def get_idx():
    print("Choose index to select variable(s) to scale/plot ,  example -hist 4   or -hist 4,5,12,11")
    for k, v in sorted(var_mapping.items(), key=lambda x : x[0]):
            print(" ##  hist {}  : {}".format(v, k))


def main(hist, idx, ch, isblinded):
    for_combine = False
    if hist >= 34 and hist<40 :
        for_combine = True
        hist = 32
    print('hist = ', hist, " : ", var_mapping[int(hist)])
    bashCommandsA =  [ 
                    "echo running {} ".format(var_mapping[int(hist)]),
                    "bash remove_existingfiles.sh scaled_files {}".format(var_mapping[int(hist)]),
                    "python3 postAnalyzer.py --hist {} --idx {} --blinded {}".format(hist, idx, isblinded),
                     ]
    bashCommandsB =  [      
                    "bash remove_existingfiles.sh agg_file {}".format(var_mapping[int(hist)]),
                    "hadd f_{i}_initial.root sample/*{i}.root".format(i=var_mapping[int(hist)]),
                    #"hadd f_{i}_up.root sample/*{i}_up.root".format(i=var_mapping[int(hist)]),
                    #"hadd f_{i}_down.root sample/*{i}_down.root".format(i=var_mapping[int(hist)]),             
                    "python3 get_zll.py --hist {}".format(var_mapping[int(hist)]),
                    "python3 get_small_mc.py --hist {}".format(var_mapping[int(hist)]),
                    "python3 get_jetFakes.py --hist {} --blinded {}".format(var_mapping[int(hist)], isblinded),
                    #"python3 get_jetFakes_unc.py --hist {} --blinded {}".format(var_mapping[int(hist)], isblinded),
                     ]
    bashCommandsC =  [ 
                    "bash remove_existingfiles.sh agg_file {}".format(ch),
                    "hadd f_{}_initial.root sample/*tot_TMass_full.root".format(ch),
                    #"hadd f_{}_up.root sample/*tot_TMass_full_up.root".format(ch),
                    #"hadd f_{}_down.root sample/*tot_TMass_full_down.root".format(ch),                    
                    "python3 get_zll.py --hist {}".format(ch),
                    "python3 get_small_mc.py --hist {}".format(ch),
                    "python3 get_jetFakes.py --hist {} --blinded {}".format(ch, isblinded),
                    #"python3 get_jetFakes_unc.py --hist {} --blinded {}".format(ch, isblinded),
                    "python3 gather_hist_v3.py",
                    "python3 gather_hist_v4.py"
                     ]
    bashCommands = []
    if  for_combine == True :
        bashCommands = bashCommandsA + bashCommandsC

    else:
        bashCommands = bashCommandsA + bashCommandsB
    command = ' ; '.join(bashCommands)
    start = time.time()
    print("executing  |  " + command)

    # try:
    #     process = subprocess.Popen(command , shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #     while process.poll() is None:
    #         print(process.stdout.readline())
    #     process.communicate()
    # except KeyboardInterrupt:
    #     print("Got Keyboard interrupt"    )
    for cmd in bashCommands:
        subprocess.check_call( cmd , shell=True)
    
    end = time.time()
    print( 'Time for execution = ', round((end - start)/60, 1) ,'minutes' )
    print("Total time for executing = {} minutes ".format((end-start)//60))
    print("....Done.....")

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    get_idx()
    print(" ")
    print("Usage : python main.py --hist 1 --idx 9 --ch etau --y 2017")
    parser.add_argument("-hist",
                    help="index of histogram to be plotted,  example -h 4   or -idx 4,5,12,11 ")       
    parser.add_argument("-idx",
                    help="index of selection to be plotted,  example -idx 9, default is 9 ")    
                    
    parser.add_argument("-ch",
                    help="channel name, etau, mutau, tautau ")
    parser.add_argument("-y",
                    help="year 2016,2017,2018")
    parser.add_argument("--blinded",
                    help="is this for blinded case using 1/5th data,   Default= 0",
                    choices=('0', '1'),
                    default='0'
                    )
    args =  parser.parse_args()
    idx = ''
    if args.y is None:
        print("Specify year, 2016,2017,2018")
        exit()     
    if args.ch is None:
        print("Specify channel name etau, mutau or tautau")
        exit()    
    
    if args.idx is None:
        if args.y=='2018':
            idx = '9b'
        else:
            idx = '9'
        print("No index passed, Specify selection index " +idx)
    else:
        idx = args.idx

    if args.hist is None:
        print("Specify histogram index ")
        exit()
    channel = args.ch
    #print int(args.hist), idx, channel, args.blinded
    for hist in args.hist.split(','):
        main(int(hist), idx, channel, args.blinded)
