from HTTutilities.Jet2TauFakes.Utilities import Leaf, Node, fill, FakeFactor, replace_nodes
import ROOT
import os

#Meta-data
version='20160914'
tag='v0.2.1'
channels=["tt"]
categories = ['incl','_0jet','_1jet','_1jetZ050','_1jetZ50100','_1jetZ100','_2jet','_2jetVBF','_anyb']

for channel in channels:
    for category in categories:

        print 'Fake factor input file for channel {0} and category {1}'.format(channel,category)
        
        # Individual fake factors
        ff_qcd_os = FakeFactor(vars=['tau_pt', 'tau_decay', 'njets', 'mvis', 'mt', 'mu_iso'])
        # Combined fake factor
        ff_comb   = FakeFactor(vars=['tau_pt', 'tau_decay', 'njets', 'mvis', 'mt', 'mu_iso'])
        
        
        home = os.getenv('HOME')
        
        ###########################################################################################################
        ### QCD fake factors
        
        qcd_os = Node(
            name='ff_qcd_os',
            formula='{mviscorr_qcd}*{ff_raw_qcd}', 
            leaves=[
                Leaf(
                    name='ff_raw_qcd',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_QCD_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='FakeFactors_Data_QCDSS_3D_SS_Iso_Medium_SS_InvertIso_Medium_tau_pt_vs_decayMode',
                    vars=['tau_pt','tau_decay','njets']
                ),
                Leaf(
                    name='mviscorr_qcd',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/Correction_Data_QCD_MVis.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='QCD_SS_MuMedium_Data_FFSSMuMediumData_mvis_correction',
                    vars=['mvis']
                )
            ]
        )
        qcd_os_up = replace_nodes(
            qcd_os, 
            {'ff_qcd_os':
             Node(
                 name='ff_qcd_os_up',
                 formula='(1.+{sys_qcd_up})*{ff_qcd_os}',
                 leaves=[
                     Leaf(
                         name='sys_qcd_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_QCD_MVis_Iso_SS2OS_up',
                         vars=['mvis']
                     ),
                     qcd_os.find('ff_qcd_os')
                 ]
             )
            }
        )
        qcd_os_down = replace_nodes(
            qcd_os, 
            {'ff_qcd_os':
             Node(
                 name='ff_qcd_os_down',
                 formula='max(0.,1.-{sys_qcd_down})*{ff_qcd_os}',
                 leaves=[
                     Leaf(
                         name='sys_qcd_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_QCD_MVis_Iso_SS2OS_down',
                         vars=['mvis']
                     ),
                     qcd_os.find('ff_qcd_os')
                 ]
             )
            }
        )
        qcd_os_up_stat = replace_nodes(
            qcd_os, 
            {'ff_qcd_os':
             Node(
                 name='ff_qcd_os_up_stat',
                 formula='(1.+{stat_qcd_up})*{ff_qcd_os}',
                 leaves=[
                     Leaf(
                         name='stat_qcd_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_QCD_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_QCDSS_3D_SS_Iso_Medium_SS_InvertIso_Medium_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     qcd_os.find('ff_qcd_os')
                 ]
             )
            }
        )
        
        qcd_os_down_stat = replace_nodes(
            qcd_os, 
            {'ff_qcd_os':
             Node(
                 name='ff_qcd_os_down_stat',
                 formula='max(0.,1.-{stat_qcd_down})*{ff_qcd_os}',
                 leaves=[
                     Leaf(
                         name='stat_qcd_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_QCD_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_QCDSS_3D_SS_Iso_Medium_SS_InvertIso_Medium_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     qcd_os.find('ff_qcd_os')
                 ]
             )
            }
        )        
        
        ###########################################################################################################
        ### Combined fake factors
        comb = Node(
            name='ff_comb',
            formula='({frac_tt}+{frac_w}+{frac_dy}+{frac_qcd})*{ff_qcd_os}',
            leaves=[
                # Individual fake factors
                qcd_os,
                # Fractions
                Leaf(
                    name='frac_qcd',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/frac_qcd.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='h_w_2d',
                    vars=['mt','tau_decay']
                ),
                Leaf(
                    name='frac_w',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/frac_wjets.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='h_w_2d',
                    vars=['mt','tau_decay']
                ),
                Leaf(
                    name='frac_dy',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/frac_dy.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='h_w_2d',
                    vars=['mt','tau_decay']
                ),
                Leaf(
                    name='frac_tt',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/frac_tt.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='h_w_2d',
                    vars=['mt','tau_decay']
                ),
            ]
        )
        
        comb_qcd_up = replace_nodes(
            comb,
            {'ff_qcd_os':
             Node(
                 name='ff_qcd_os_up',
                 formula='(1.+{sys_qcd_up})*{ff_qcd_os}',
                 leaves=[
                     Leaf(
                         name='sys_qcd_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_QCD_MVis_Iso_SS2OS_up',
                         vars=['mvis']
                     ),
                     qcd_os.find('ff_qcd_os')
                 ]
             )
            }
        )
        comb_qcd_down = replace_nodes(
            comb,
            {'ff_qcd_os':
            Node(
                name='ff_qcd_os_down',
                formula='max(0.,1.-{sys_qcd_down})*{ff_qcd_os}',
                leaves=[
                    Leaf(
                        name='sys_qcd_down',
                        file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                        object='uncertainties_QCD_MVis_Iso_SS2OS_down',
                        vars=['mvis']
                    ),
                    qcd_os.find('ff_qcd_os')
                ]
            )
            }
        )
        comb_qcd_up_stat = replace_nodes(
            comb, 
            {'ff_qcd_os':
             Node(
                 name='ff_qcd_os_up_stat',
                 formula='(1.+{stat_qcd_up})*{ff_qcd_os}',
                 leaves=[
                     Leaf(
                         name='stat_qcd_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_QCD_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_QCDSS_3D_SS_Iso_Medium_SS_InvertIso_Medium_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     qcd_os.find('ff_qcd_os')
                 ]
             )
            }
        )
        comb_qcd_down_stat = replace_nodes(
            comb, 
            {'ff_qcd_os':
             Node(
                 name='ff_qcd_os_down_stat',
                 formula='max(0.,1.-{stat_qcd_down})*{ff_qcd_os}',
                 leaves=[
                     Leaf(
                         name='stat_qcd_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_QCD_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_QCDSS_3D_SS_Iso_Medium_SS_InvertIso_Medium_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     qcd_os.find('ff_qcd_os')
                 ]
             )
            }
        )        
        fill(ff_qcd_os, qcd_os)
        fill(ff_qcd_os, qcd_os_up,   sys='ff_qcd_syst_up')
        fill(ff_qcd_os, qcd_os_down, sys='ff_qcd_syst_down')
        fill(ff_qcd_os, qcd_os_up_stat,   sys='ff_qcd_stat_up')
        fill(ff_qcd_os, qcd_os_down_stat, sys='ff_qcd_stat_down')
        fill(ff_comb  , comb)
        fill(ff_comb, comb_qcd_up,   sys='ff_qcd_syst_up')
        fill(ff_comb, comb_qcd_down,   sys='ff_qcd_syst_down')
        fill(ff_comb, comb_qcd_up_stat,   sys='ff_qcd_stat_up')
        fill(ff_comb, comb_qcd_down_stat,   sys='ff_qcd_stat_down')
        
        
        
        file = ROOT.TFile.Open("{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/fakeFactors_{VERSION}.root".format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category), "recreate")
        # Write meta-data
        version_ts = ROOT.TString(version)
        tag_ts     = ROOT.TString(tag)
        file.WriteObject(version_ts , "version")
        file.WriteObject(tag_ts     , "tag")
        # Write fake factors
        file.WriteObject(ff_qcd_os.fakefactor  , "ff_qcd_os")
        file.WriteObject(ff_comb.fakefactor    , "ff_comb")
        #
        file.Close()
    
