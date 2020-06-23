import os
sampleListDict = {
    'DY1JetsToLL' : ['DY1JetsToLL', 'ZTT', '0'], 
    'DY2JetsToLL' : ['DY2JetsToLL', 'ZTT', '0'], 
    'DY3JetsToLL' : ['DY3JetsToLL', 'ZTT', '0'], 
    'DY4JetsToLL' : ['DY4JetsToLL', 'ZTT', '0'], 
    'DYJetsToLL'  : ['DYJetsToLL' , 'ZTT', '0'], 
    'DYJetsToLL_M10to50' : ['DYJetsToLL_M10to50', 'ZTT', '0'], 
    
    'DY1JetsToLL_stitch' : ['ZTT1jet', 'ZTTjet', '0'],
    'DY2JetsToLL_stitch' : ['ZTT2jet', 'ZTTjet', '0'],
    'DY3JetsToLL_stitch' : ['ZTT3jet', 'ZTTjet', '0'],
    'DY4JetsToLL_stitch' : ['ZTT4jet', 'ZTTjet', '0'],
    'DYJetsToLL_stitch'  : ['ZTTjet_inc' , 'ZTTjet', '0'],

    'EWKWMinus2Jets_WToLNu' : ['EWKWMinus2Jets', 'EWKWMinus', '0'], 
    'EWKWPlus2Jets_WToLNu' : ['EWKWPlus2Jets', 'EWKWPlus', '0'], 
    'EWKZ2Jets_ZToLL' : ['EWKZ2Jets_ZToLL', 'EWKZ2Jets', '0'], 
    'EWKZ2Jets_ZToNuNu' : ['EWKZ2Jets_ZToNuNu', 'EWKZ2Jets', '0'], 
    'GluGluHToTauTau' : ['GluGluHToTauTau', 'GluGluH', '0'], 
    'GluGluHToWWTo2L2Nu' : ['GluGluHToWWTo2L2Nu', 'GluGluH', '0'], 
    'GluGluZH_HToWW_' : ['GluGluZH_HToWW', 'GluGluZH', '0'], 
    'HWminusJ_HToWW' : ['HWminusJ_HToWW', 'HWminusJ', '0'], 
    'HWplusJ_HToWW' : ['HWplusJ_HToWW', 'HWplusJ', '0'], 
    'ZZTo4L' : ['ZZTo4L', 'VV', '0'], 
    'HZJ_HToWW' : ['HZJ_HToWW', 'HZJ', '0'], 
    'ST_t-channel_antitop' : ['ST_t-channel_antitop', 'ST_t', '0'], 
    'ST_t-channel_top' : ['ST_t-channel_top', 'ST_t', '0'], 
    'ST_tW_antitop' : ['ST_tW_antitop', 'ST_t', '0'], 
    'ST_tW_top' : ['ST_tW_top', 'ST_t', '0'], 
    'SingleMuonA' : ['data_obs', 'data_obs', '0'], 
    'SingleMuonB' : ['data_obs', 'data_obs', '0'], 
    'SingleMuonC' : ['data_obs', 'data_obs', '0'], 
    'SingleMuonD_PromptReco' : ['data_obs', 'data_obs', '0'], 
    'TTTo2L2Nu_powheg' : ['TTTo2L2Nu', 'TT', '0'], 
    'TTToHadronic_powheg' : ['TTToHadronic', 'TT', '0'], 
    'TTToSemiLeptonic_powheg' : ['TTToSemiLeptonic', 'TT', '0'], 
    'VBFHToTauTau' : ['VBFHToTauTau', 'VBFH', '0'], 
    'VBFHToWWTo2L2Nu' : ['VBFHToWWTo2L2Nu', 'VBFH', '0'], 
    'VVTo2L2Nu' : ['VVTo2L2Nu', 'VV', '0'], 
    'W1JetsToLNu' : ['W1JetsToLNu', 'WJets', '0'], 
    'W2JetsToLNu' : ['W2JetsToLNu', 'WJets', '0'], 
    'W3JetsToLNu' : ['W3JetsToLNu', 'WJets', '0'], 
    'W4JetsToLNu' : ['W4JetsToLNu', 'WJets', '0'],
    'WJetsToLNu' : ['WJetsToLNu', 'WJets', '0'],
    
    'W1JetsToLNu_stitch' : ['W1Jet', 'WJets_jets', '0'],
    'W2JetsToLNu_stitch' : ['W2Jet', 'WJets_jets', '0'],
    'W3JetsToLNu_stitch' : ['W3Jet', 'WJets_jets', '0'],
    'W4JetsToLNu_stitch' : ['W4Jet', 'WJets_jets', '0'],
    'WJetsToLNu_stitch'  : ['WJets_inc' , 'WJets_jets', '0'],

    'WGToLNuG' : ['WGToLNuG', 'WGToLNuG', '0'], 
    'WWTo1L1Nu2Q' : ['WWTo1L1Nu2Q', 'VV', '0'], 
    'WWToLNuQQ' : ['WWToLNuQQ', 'VV', '0'], 
    'WZTo3LNu' : ['WZTo3LNu', 'VV', '0'], 
    'WminusHToTauTau' : ['WminusHToTauTau', 'WminusH', '0'], 
    'WplusHToTauTau' : ['WplusHToTauTau', 'WplusH', '0'], 
    'ZHToTauTau' : ['ZHToTauTau', 'ZH', '0'], 
    'ZZTo2L2Q' : ['ZZTo2L2Q', 'VV', '0'], 
    'ggZH_HToTauTau_ZToLL' : ['ggZH_HToTauTau_ZToLL', 'ggZH', '0'], 
    'ggZH_HToTauTau_ZToNuNu' : ['ggZH_HToTauTau_ZToNuNu', 'ggZH', '0'], 
    'ggZH_HToTauTau_ZToQQ' : ['ggZH_HToTauTau_ZToQQ', 'ggZH', '0'], 
    'ttHToNonbb' : ['ttHToNonbb', 'ttH', '0'], 
    'WWW' :  ['WWW' , 'VVV' ,'0' ],
    'WWZ' :  ['WWZ' , 'VVV' ,'0' ],
    'WZZ' :  ['WZZ' , 'VVV' ,'0' ],
    'ZZZ' :  ['ZZZ' , 'VVV' ,'0' ],
    'ZJetsToNuNu_HT100-200' :   ['ZJetsToNuNu_HT100-200'   , 'ZJetsToNuNu' ,'0'],
    'ZJetsToNuNu_HT1200-2500' : ['ZJetsToNuNu_HT1200-2500' , 'ZJetsToNuNu' ,'0'],
    'ZJetsToNuNu_HT200-400' :   ['ZJetsToNuNu_HT200-400'   , 'ZJetsToNuNu' ,'0'],
    'ZJetsToNuNu_HT2500-Inf' :  ['ZJetsToNuNu_HT2500-Inf'  , 'ZJetsToNuNu' ,'0'],
    'ZJetsToNuNu_HT400-600' :   ['ZJetsToNuNu_HT400-600'   , 'ZJetsToNuNu' ,'0'],
    'ZJetsToNuNu_HT600-800' :   ['ZJetsToNuNu_HT600-800'   , 'ZJetsToNuNu' ,'0'],
    'ZJetsToNuNu_HT800-1200' :  ['ZJetsToNuNu_HT800-1200'  , 'ZJetsToNuNu' ,'0' ],
    
}

filelist=os.listdir("files_initial")
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
#sh do_plots_mutau.sh 
echo "*************** plots made ***************" 
"""
)

outFile.close()

