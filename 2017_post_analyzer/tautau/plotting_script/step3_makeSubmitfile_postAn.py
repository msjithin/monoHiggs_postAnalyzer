import os
sampleListDict = {
    'DY1JetsToLL_M-50_TuneCP5' : ['DY1JetsToLL', 'ZTT', '0'], 
    'DY2JetsToLL_M-50_TuneCP5' : ['DY2JetsToLL', 'ZTT', '0'], 
    'DY3JetsToLL_M-50_TuneCP5_ext1' : ['DY3JetsToLL', 'ZTT', '0'], 
    'DY3JetsToLL_M-50_TuneCP5_v1' : ['DY3JetsToLL', 'ZTT', '0'], 
    'DY4JetsToLL_M-50_TuneCP5' : ['DY4JetsToLL', 'ZTT', '0'], 
    'DYJetsToLL_M-10to50_TuneCP5' : ['DYJetsToLL', 'ZTT', '0'], 
    'DYJetsToLL_M-50_TuneCP5_ext1_v1' : ['DYJetsToLL', 'ZTT', '0'], 
    'DYJetsToLL_M-50_TuneCP5_v1' : ['DYJetsToLL', 'ZTT', '0'], 
    
    'DY1JetsToLL_M-50_TuneCP5_stitch' : ['ZTT1jet', 'ZTTjet', '0'],
    'DY2JetsToLL_M-50_TuneCP5_stitch' : ['ZTT2jet', 'ZTTjet', '0'],
    'DY3JetsToLL_M-50_TuneCP5_ext1_stitch' : ['ZTT3jet', 'ZTTjet', '0'],
    'DY3JetsToLL_M-50_TuneCP5_v1_stitch' : ['ZTT3jet', 'ZTTjet', '0'],
    'DY4JetsToLL_M-50_TuneCP5_stitch' : ['ZTT4jet', 'ZTTjet', '0'],
    'DYJetsToLL_M-50_TuneCP5_ext1_v1_stitch'  : ['ZTTjet_inc' , 'ZTTjet', '0'],
    'DYJetsToLL_M-50_TuneCP5_v1_stitch'  : ['ZTTjet_inc' , 'ZTTjet', '0'],

    'EWKWMinus2Jets_WToLNu_M-50_TuneCP5' : ['EWKWMinus2Jets', 'EWKWMinus', '0'], 
    'EWKWPlus2Jets_WToLNu_M-50_TuneCP5' : ['EWKWPlus2Jets', 'EWKWPlus', '0'], 
    'EWKZ2Jets_ZToLL_M-50_TuneCP5' : ['EWKZ2Jets_ZToLL', 'EWKZ2Jets', '0'], 
    'EWKZ2Jets_ZToNuNu_TuneCP5' : ['EWKZ2Jets_ZToNuNu', 'EWKZ2Jets', '0'], 
    'GluGluHToTauTau_M125' : ['GluGluHToTauTau', 'GluGluH', '0'], 
    'GluGluHToWWTo2L2Nu_M125' : ['GluGluHToWWTo2L2Nu', 'GluGluH', '0'], 
    'ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5' : ['ST_t-channel_antitop', 'ST_t', '0'], 
    'ST_t-channel_top_4f_inclusiveDecays_TuneCP5' : ['ST_t-channel_top', 'ST_t', '0'], 
    'ST_tW_antitop_5f_inclusiveDecays_TuneCP5' : ['ST_tW_antitop', 'ST_t', '0'], 
    'ST_tW_top_5f_inclusiveDecays_TuneCP5' : ['ST_tW_top', 'ST_t', '0'], 
    'Tau_EraB' : ['data_obs', 'data_obs', '0'], 
    'Tau_EraC' : ['data_obs', 'data_obs', '0'], 
    'Tau_EraD' : ['data_obs', 'data_obs', '0'], 
    'Tau_EraE' : ['data_obs', 'data_obs', '0'], 
    'Tau_EraF' : ['data_obs', 'data_obs', '0'], 
    'TTTo2L2Nu_TuneCP5' : ['TTTo2L2Nu', 'TT', '0'], 
    'TTToHadronic_TuneCP5' : ['TTToHadronic', 'TT', '0'], 
    'TTToSemiLeptonic_TuneCP5' : ['TTToSemiLeptonic', 'TT', '0'], 
    'VBFHToTauTau_M125' : ['VBFHToTauTau', 'VBFH', '0'], 
    'VBFHToWWTo2L2Nu_M125' : ['VBFHToWWTo2L2Nu', 'VBFH', '0'], 
    'VVTo2L2Nu' : ['VVTo2L2Nu', 'VV', '0'], 
    'W1JetsToLNu_TuneCP5' : ['W1JetsToLNu', 'WJets', '0'], 
    'W2JetsToLNu_TuneCP5' : ['W2JetsToLNu', 'WJets', '0'], 
    'W3JetsToLNu_TuneCP5' : ['W3JetsToLNu', 'WJets', '0'], 
    'W4JetsToLNu_TuneCP5' : ['W4JetsToLNu', 'WJets', '0'], 
    'WJetsToLNu_TuneCP5' : ['WJetsToLNu', 'WJets', '0'], 
    
    'W1JetsToLNu_TuneCP5_stitch' : ['W1Jet', 'WJets_jets', '0'],
    'W2JetsToLNu_TuneCP5_stitch' : ['W2Jet', 'WJets_jets', '0'],
    'W3JetsToLNu_TuneCP5_stitch' : ['W3Jet', 'WJets_jets', '0'],
    'W4JetsToLNu_TuneCP5_stitch' : ['W4Jet', 'WJets_jets', '0'],
    'WJetsToLNu_TuneCP5_stitch'  : ['WJets_inc' , 'WJets_jets', '0'],
    
    'WWTo1L1Nu2Q' : ['WWTo1L1Nu2Q', 'VV', '0'], 
    'WWToLNuQQ_NNPDF31_TuneCP5' : ['WWToLNuQQ', 'VV', '0'], 
    'WWW_4F_TuneCP5' : ['WWW', 'VVV', '0'], 
    'WWZ_4F_TuneCP5' : ['WWZ', 'VVV', '0'], 
    'WW_TuneCP5' : ['WW', 'VV', '0'], 
    'WZTo3LNu_TuneCP5' : ['WZTo3LNu', 'VV', '0'], 
    'WZZ_TuneCP5' : ['WZZ', 'VVV', '0'], 
    'WZ_TuneCP5' : ['WZ', 'VV', '0'], 
    'WminusHToTauTau_M125' : ['WminusHToTauTau', 'WminusH', '0'], 
    'WplusHToTauTau_M125' : ['WplusHToTauTau', 'WplusH', '0'], 
    'ZHToTauTau_M125' : ['ZHToTauTau', 'ZH', '0'], 
    'ZJetsToNuNu_HT-100To200' : ['ZJetsToNuNu_HT100-200', 'ZJetsToNuNu', '0'], 
    'ZJetsToNuNu_HT-1200To2500' : ['ZJetsToNuNu_HT1200-2500', 'ZJetsToNuNu', '0'], 
    'ZJetsToNuNu_HT-200To400' : ['ZJetsToNuNu_HT200-400', 'ZJetsToNuNu', '0'], 
    'ZJetsToNuNu_HT-2500ToInf' : ['ZJetsToNuNu_HT2500-Inf', 'ZJetsToNuNu', '0'], 
    'ZJetsToNuNu_HT-400To600' : ['ZJetsToNuNu_HT400-600', 'ZJetsToNuNu', '0'], 
    'ZJetsToNuNu_HT-600To800' : ['ZJetsToNuNu_HT600-800', 'ZJetsToNuNu', '0'], 
    'ZJetsToNuNu_HT-800To1200' : ['ZJetsToNuNu_HT800-1200', 'ZJetsToNuNu', '0'], 
    'ZZTo2L2Q' : ['ZZTo2L2Q', 'VV', '0'], 
    'ZZTo4L_TuneCP5' : ['ZZTo4L', 'VV', '0'], 
    'ZZZ_TuneCP5' : ['ZZZ', 'VVV', '0'], 
    'ZZ_TuneCP5' : ['ZZ', 'VV', '0'], 
    'ggZH_HToTauTau_ZToLL_M125' : ['ggZH_HToTauTau_ZToLL', 'ggZH', '0'], 
    'ggZH_HToTauTau_ZToNuNu_M125' : ['ggZH_HToTauTau_ZToNuNu', 'ggZH', '0'], 
    'ttHToTauTau_M125_TuneCP5' : ['ttHToTauTau', 'ttH', '0'], 
    
}

filelist=os.listdir("../files_initial")
filelist = [item for item in filelist if '.root' in item]

samplelist=[]
samplelist = list(map(lambda x: 
                      x.replace('.root','')
                      .replace('_final','')
                      , filelist))

#for n,g in zip(filelist, samplelist):
#    print(n+'\t\t\t'+g)
outFile = open("_postAnalyzer_mutau.sh", "w")
outFile.write(
"""
#!/bin/bash 
set -e 
if [ -f "f_mutau_initial.root" ]; then
    echo "deleting existing f_mutau_initial.root file ....."
    rm f_mutau_initial.root
fi
if [ "$(ls -A files_nominal)" ]; then
    echo "deleting existing files in directory files_nominal ....."
    rm files_nominal/*.root
fi


./Make.sh _postAnalyzer_mutau.C 

"""
)
for j in range(len(filelist)) :
    tmp = samplelist[j]
    samplename=sampleListDict[tmp]
    outoutFilename=filelist[j]
    outFile.write('./_postAnalyzer_mutau.exe ../files_initial/{} files_nominal/{} {} {} {} \n'
                  .format(filelist[j], outoutFilename, samplename[0], samplename[1], samplename[2]))
    st_tmp=tmp+'_stitch'
    if st_tmp in sampleListDict :
        samplename=sampleListDict[st_tmp]
        outoutFilename=samplelist[j]+'_stitch_final.root'
        outFile.write('./_postAnalyzer_mutau.exe ../files_initial/{} files_nominal/{} {} {} {} \n'
                      .format(filelist[j], outoutFilename, samplename[0], samplename[1], samplename[2]))

outFile.write(
"""
hadd -f f_mutau_initial.root files_nominal/*.root 
echo "*************** root file made ***************" 
sh do_make_plots.sh
echo "*************** plots made ***************" 
"""
)

outFile.close()

print("""
output written to _postAnalyzer_mutau.sh


""")
