#!/usr/bin/env python
import ROOT
from os import path, remove, getcwd

mc_samples = ['jetFakes', 'ZTTjet', 'ZLLjet', 'TT' , 'otherMC' , 'STT', 'VVT']
# signal_samples=['MH3_200_MH4_100', 'MH3_200_MH4_150', 'MH3_300_MH4_100', 'MH3_300_MH4_150', 'MH3_400_MH4_100', 'MH3_400_MH4_150', 'MH3_400_MH4_200', 'MH3_400_MH4_250', 'MH3_500_MH4_150', 'MH3_500_MH4_200', 'MH3_500_MH4_250', 'MH3_500_MH4_300', 'MH3_600_MH4_100', 'MH3_600_MH4_150', 'MH3_600_MH4_200', 'MH3_600_MH4_250', 'MH3_600_MH4_300', 'MH3_600_MH4_350', 'MH3_600_MH4_400', 'MH3_600_MH4_500', 'MH3_700_MH4_250', 'MH3_700_MH4_300', 'MH3_700_MH4_350', 'MH3_700_MH4_400', 'MH3_800_MH4_250', 'MH3_800_MH4_300', 'MH3_800_MH4_350', 'MH3_800_MH4_500', 'MH3_900_MH4_300', 'MH3_900_MH4_350', 'MH3_900_MH4_400', 'MH3_900_MH4_500']
signal_samples=['2HDMa_bb_sinp_0p1_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_bb_sinp_0p2_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_bb_sinp_0p35_tanb_0p5_mXd_10_MH3_600_MH4_150', '2HDMa_bb_sinp_0p35_tanb_0p5_mXd_10_MH3_600_MH4_250', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_1000_MH4_150', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_1000_MH4_350', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_1200_MH4_150', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_1200_MH4_250', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_1200_MH4_350', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_1600_MH4_150', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_1600_MH4_250', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_1600_MH4_350', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_200_MH4_150', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_400_MH4_150', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_400_MH4_250', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_150', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_250', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_350', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_800_MH4_250', '2HDMa_bb_sinp_0p35_tanb_1p0_mXd_10_MH3_800_MH4_350', '2HDMa_bb_sinp_0p35_tanb_1p5_mXd_10_MH3_600_MH4_150', '2HDMa_bb_sinp_0p35_tanb_1p5_mXd_10_MH3_600_MH4_250', '2HDMa_bb_sinp_0p35_tanb_20p0_mXd_10_MH3_600_MH4_150', '2HDMa_bb_sinp_0p35_tanb_20p0_mXd_10_MH3_600_MH4_250', '2HDMa_bb_sinp_0p35_tanb_2p0_mXd_10_MH3_600_MH4_150', '2HDMa_bb_sinp_0p35_tanb_2p0_mXd_10_MH3_600_MH4_250', '2HDMa_bb_sinp_0p35_tanb_4p0_mXd_10_MH3_600_MH4_150', '2HDMa_bb_sinp_0p35_tanb_4p0_mXd_10_MH3_600_MH4_250', '2HDMa_bb_sinp_0p35_tanb_50p0_mXd_10_MH3_600_MH4_150', '2HDMa_bb_sinp_0p35_tanb_50p0_mXd_10_MH3_600_MH4_250', '2HDMa_bb_sinp_0p35_tanb_8p0_mXd_10_MH3_600_MH4_150', '2HDMa_bb_sinp_0p35_tanb_8p0_mXd_10_MH3_600_MH4_250', '2HDMa_bb_sinp_0p3_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_bb_sinp_0p4_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_bb_sinp_0p6_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_bb_sinp_0p7_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_bb_sinp_0p8_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_bb_sinp_0p9_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p1_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p2_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p35_tanb_0p5_mXd_10_MH3_600_MH4_150', '2HDMa_gg_sinp_0p35_tanb_0p5_mXd_10_MH3_600_MH4_250', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_1000_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_1000_MH4_250', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_1000_MH4_350', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_1200_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_1200_MH4_250', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_1200_MH4_350', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_1600_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_1600_MH4_250', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_1600_MH4_350', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_800_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p5_mXd_10_MH3_600_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p5_mXd_10_MH3_600_MH4_250', '2HDMa_gg_sinp_0p35_tanb_20p0_mXd_10_MH3_600_MH4_150', '2HDMa_gg_sinp_0p35_tanb_20p0_mXd_10_MH3_600_MH4_250', '2HDMa_gg_sinp_0p35_tanb_2p0_mXd_10_MH3_600_MH4_150', '2HDMa_gg_sinp_0p35_tanb_2p0_mXd_10_MH3_600_MH4_250', '2HDMa_gg_sinp_0p35_tanb_4p0_mXd_10_MH3_600_MH4_150', '2HDMa_gg_sinp_0p35_tanb_4p0_mXd_10_MH3_600_MH4_250', '2HDMa_gg_sinp_0p35_tanb_50p0_mXd_10_MH3_600_MH4_150', '2HDMa_gg_sinp_0p35_tanb_50p0_mXd_10_MH3_600_MH4_250', '2HDMa_gg_sinp_0p35_tanb_8p0_mXd_10_MH3_600_MH4_150', '2HDMa_gg_sinp_0p35_tanb_8p0_mXd_10_MH3_600_MH4_250', '2HDMa_gg_sinp_0p3_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p4_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p5_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p6_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p7_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p8_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p9_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_200_MH4_100', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_200_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_300_MH4_100', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_300_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_400_MH4_100', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_400_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_400_MH4_200', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_400_MH4_250', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_500_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_500_MH4_200', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_500_MH4_250', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_500_MH4_300', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_100', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_150', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_200', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_250', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_300', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_350', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_400', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_500', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_700_MH4_250', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_700_MH4_300', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_700_MH4_350', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_700_MH4_400', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_800_MH4_250', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_800_MH4_300', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_800_MH4_350', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_800_MH4_500', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_900_MH4_300', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_900_MH4_350', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_900_MH4_400', '2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_900_MH4_500']

zprime_samples=['MZp_1000_MChi_100', 'MZp_1000_MChi_1', 'MZp_1000_MChi_200', 'MZp_1000_MChi_400', 'MZp_1000_MChi_600', 'MZp_1000_MChi_800', 'MZp_100_MChi_1', 'MZp_100_MChi_50', 'MZp_1500_MChi_100', 'MZp_1500_MChi_1', 'MZp_1500_MChi_200', 'MZp_1500_MChi_400', 'MZp_1500_MChi_600', 'MZp_1500_MChi_800', 'MZp_2000_MChi_100', 'MZp_2000_MChi_1', 'MZp_2000_MChi_200', 'MZp_2000_MChi_400', 'MZp_2000_MChi_600', 'MZp_2000_MChi_800', 'MZp_200_MChi_100', 'MZp_200_MChi_150', 'MZp_200_MChi_1', 'MZp_200_MChi_50', 'MZp_2500_MChi_100', 'MZp_2500_MChi_1', 'MZp_2500_MChi_200', 'MZp_2500_MChi_400', 'MZp_2500_MChi_600', 'MZp_2500_MChi_800', 'MZp_3000_MChi_100', 'MZp_3000_MChi_1', 'MZp_3000_MChi_200', 'MZp_300_MChi_150', 'MZp_3500_MChi_100', 'MZp_3500_MChi_1', 'MZp_350_MChi_50', 'MZp_500_MChi_100', 'MZp_500_MChi_1', 'MZp_500_MChi_200', 'MZp_500_MChi_400', 'MZp_650_MChi_50', 'MZp_800_MChi_50']


        
def checkHistogram(f, histogram):
    isthere=  f.GetListOfKeys().Contains(histogram)
    #print(isthere)
    return isthere


def GetKeyNames( self, dir = "" ):
    self.cd(dir)
    return [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
ROOT.TFile.GetKeyNames = GetKeyNames

def getLabel(name, idx):
    label = name
    label = label.replace('tot_TMass_full_'+idx, '')
    label = label.replace('_down1', 'Down')
    label = label.replace('_up1', 'Up')
    label = label.replace('_down', 'Down')
    label = label.replace('_up', 'Up')
    label = label.replace('_Down', 'Down')
    label = label.replace('_Up', 'Up')
    if '_fr' in label:
        label = label.replace('_fr', '')
    return label

def getHistList(fname , channelName, idx):
    #channelName = 'etau'
    outFile =  ROOT.TFile(channelName+'.root',"UPDATE")
    inFile = ROOT.TFile(fname,"READ")
    keyList = inFile.GetKeyNames()
    #print "\nKeys in file:", keyList
    tmpList= []
    channel = channelName
    outFile.cd()
    if not outFile.GetDirectory(channel):
        outFile.mkdir(channel)
    outFile.cd(channel)

    for tdir in sorted(keyList):
        if 'initial' in fname and tdir!='tot_TMass_full_'+idx:
            continue
        elif ('up' in fname or 'down' in fname) and 'tot_TMass_full_'+idx+'_' not in tdir:
            continue
        #if 'tot_TMass_full' not in tdir : continue
        #if '_9' not in tdir: continue
        if "_dyll" in tdir: continue
        if tdir == 'tot_TMass_full_'+idx:
            data = inFile.Get(tdir+'/data_obs_'+tdir)
            data.SetName('data_obs')
            outFile.cd(channel)
            data.Write()

            if inFile.Get(tdir+'/blinded_data_obs_'+tdir):
                full_lumi_data = inFile.Get(tdir+'/blinded_data_obs_'+tdir)
                full_lumi_data.SetName('blinded_data_obs')
                outFile.cd(channel)
                full_lumi_data.Write()
            else:
                print("You didnt have full luminosity data")
        for signal in signal_samples:
            if inFile.Get(tdir+'/'+signal+'_'+tdir) and '_fr' not in tdir:
                signalhist = inFile.Get(tdir+'/'+signal+'_'+tdir)
                signalhist.SetName( signal + getLabel(tdir, idx))
                signalhist.Write()
        for mc in mc_samples:
            if mc != 'jetFakes' and '_fr' in tdir:
                continue
            tmppath = tdir+'/'+mc+'_'+tdir
            try:
                tmpHist = inFile.Get(tmppath)
                tmpHist.SetName(mc + getLabel(tdir, idx))
                tmpHist.Write()
            except:
                pass
        for signal in zprime_samples:
            if inFile.Get(tdir+'/'+signal+'_'+tdir) and '_fr' not in tdir:
                signalhist = inFile.Get(tdir+'/'+signal+'_'+tdir)
                signalhist.SetName( signal + getLabel(tdir, idx))
                signalhist.Write()

    inFile.Close()
    outFile.Close()
    print("Finished gathering from " + fname)




###########################################################################################
cwd = getcwd()
#print cwd
ch = cwd.split('/')[-2]
ch = ch.replace('_blinded', '')

channelName = 'mutau'
idx = '9'
# inFile_nominal= ROOT.TFile("f_mutau_initial.root","UPDATE")
file_list = ["f_"+channelName+"_initial.root", "f_"+channelName+"_up.root", "f_"+channelName+"_down.root"]

if path.exists(channelName+".root"):
    remove(channelName+".root")
    print("The file has been deleted successfully")
else:
    print("The file does not exist!")
    
for fname in file_list:
    getHistList(fname , channelName, idx )
    
print("Done")
