from HTTutilities.Jet2TauFakes.Utilities import Leaf, Node, fill, FakeFactor, replace_nodes
import ROOT
import os

#Meta-data
version='20161023'
tag='v0.2.1'
channels=["mt","et"]
categories = ['_0jetLow','_0jetHigh','_1jetLow','_1jetHigh','_vbfLow','_vbfHigh','_2jet','_anyb']

for channel in channels:
    for category in categories:

        print 'Fake factor input file for channel {0} and category {1}'.format(channel,category)
        
        # Individual fake factors
        ff_qcd_os = FakeFactor(vars=['tau_pt', 'tau_decay', 'njets', 'mvis', 'mt', 'mu_iso'])
        ff_qcd_ss = FakeFactor(vars=['tau_pt', 'tau_decay', 'mvis', 'mu_iso'])
        ff_w      = FakeFactor(vars=['tau_pt', 'tau_decay', 'njets', 'mvis', 'mt'])
        ff_tt     = FakeFactor(vars=['tau_pt', 'tau_decay', 'njets', 'mvis', 'mt'])
        # Combined fake factor
        ff_comb   = FakeFactor(vars=['tau_pt', 'tau_decay', 'njets', 'mvis', 'mt', 'mu_iso'])
        
        
        home = os.getenv('HOME')
        
        ###########################################################################################################
        ### QCD fake factors
        
        qcd_os = Node(
            name='ff_qcd_os',
            formula='{isocorr_qcd}*{mviscorr_qcd}*{ff_raw_qcd}*{OSSS_corr_qcd}', 
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
                ),
                Leaf(
                    name='isocorr_qcd',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/Correction_Data_QCD_MuIso.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='QCD_SS_Data_FFSSMuMediumData_isomu_correction',
                    vars=['mu_iso']
                ),
                Leaf(
                    name='OSSS_corr_qcd',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/Correction_Data_QCD_OSSS.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='QCD_SS_Data_FFSSMuMediumData_OSSS_correction',
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
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_QCD_MVis_Iso_SS2OS_up',
                         vars=['mvis', 'mu_iso']
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
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_QCD_MVis_Iso_SS2OS_down',
                         vars=['mvis', 'mu_iso']
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
        ### W fake factors
        
        w = Node(
            name='ff_w',
            formula='{mtcorr_w}*{ff_raw_w}*{mviscorr_w}',
            leaves=[
                Leaf(
                    name='ff_raw_w',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_W_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='FakeFactors_Data_HighMT_3D_Iso_Medium_InvertIso_Medium_tau_pt_vs_decayMode',
                    vars=['tau_pt','tau_decay','njets']
                ),
                Leaf(
                    name='mviscorr_w',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/Correction_Data_W_MVis.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='W_OS_Data_FFOSData_mvis_correction',
                    vars=['mvis']
                ),
                Leaf(
                    name='mtcorr_w',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/Correction_MC_W_MT.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='W_OS_MC_FFOSMC_mt_correction',
                    vars=['mt']
                )
            ]
        )
        w_up = replace_nodes(
            w, 
            {'ff_w':
             Node(
                 name='ff_w_up',
                 formula='(1.+{sys_w_up})*{ff_w}',
                 leaves=[
                     Leaf(
                         name='sys_w_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_W_MVis_MT_up',
                         vars=['mvis', 'mt']
                     ),
                     w.find('ff_w')
                 ]
             )
            }
        )
        w_down = replace_nodes(
            w, 
            {'ff_w':
             Node(
                 name='ff_w_down',
                 formula='max(0.,1.-{sys_w_down})*{ff_w}',
                 leaves=[
                     Leaf(
                         name='sys_w_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_W_MVis_MT_down',
                         vars=['mvis', 'mt']
                     ),
                     w.find('ff_w')
                 ]
             )
            }
        )
        w_up_stat = replace_nodes(
            w, 
            {'ff_w':
             Node(
                 name='ff_w_up_stat',
                 formula='(1.+{stat_w_up})*{ff_w}',
                 leaves=[
                     Leaf(
                         name='stat_w_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_W_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_HighMT_3D_Iso_Medium_InvertIso_Medium_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     w.find('ff_w')
                 ]
             )
            }
        )
        w_down_stat = replace_nodes(
            w, 
            {'ff_w':
             Node(
                 name='ff_w_down_stat',
                 formula='max(0.,1.-{stat_w_down})*{ff_w}',
                 leaves=[
                     Leaf(
                         name='stat_w_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_W_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_HighMT_3D_Iso_Medium_InvertIso_Medium_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     w.find('ff_w')
                 ]
             )
            }
        )

        ###########################################################################################################
        ### TTbar fake factors
        tt = Node(
            name='ff_tt',
            formula='{ff_raw_tt}*{mviscorr_tt}',
            leaves=[
                Leaf(
                    name='ff_raw_tt',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_TT_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='FakeFactors_Data_TT_anyb_addLep_InvertIso_tau_pt_vs_decayMode',
                    vars=['tau_pt','tau_decay','njets']
                ),
                Leaf(
                    name='mviscorr_tt',
                    file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/Correction_MC_TT_MVis.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                    object='TT_OS_MC_mvis_correction',
                    vars=['mvis']
                ),
            ]
        )        
        tt_up = replace_nodes(
            tt,
            {'ff_tt':
             Node(
                 name='ff_tt_up',
                 formula='(1.+{sys_tt_up})*{ff_tt}',
                 leaves=[
                     Leaf(
                         name='sys_tt_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/uncertainties_TT.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_TT_MVis_up',
                         vars=['tau_decay','mvis']
                     ),
                     tt.find('ff_tt')
                 ]
             )
            }
        )
        tt_down = replace_nodes(
            tt,
            {'ff_tt':
             Node(
                 name='ff_tt_down',
                 formula='max(0.,1.-{sys_tt_down})*{ff_tt}',
                 leaves=[
                     Leaf(
                         name='sys_tt_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/uncertainties_TT.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_TT_MVis_down',
                         vars=['tau_decay','mvis']
                     ),
                     tt.find('ff_tt')
                 ]
             )
            }
        )
        tt_up_stat = replace_nodes(
            tt, 
            {'ff_tt':
             Node(
                 name='ff_tt_up_stat',
                 formula='(1.+{stat_tt_up})*{ff_tt}',
                 leaves=[
                     Leaf(
                         name='stat_tt_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_TT_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_TT_anyb_addLep_InvertIso_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     tt.find('ff_tt')
                 ]
             )
            }
        )
        tt_down_stat = replace_nodes(
            tt, 
            {'ff_tt':
             Node(
                 name='ff_tt_down_stat',
                 formula='max(0.,1.-{stat_tt_down})*{ff_tt}',
                 leaves=[
                     Leaf(
                         name='stat_tt_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_TT_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_TT_anyb_addLep_InvertIso_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     tt.find('ff_tt')
                 ]
             )
            }
        )
                
        ###########################################################################################################
        ### Combined fake factors
        comb = Node(
            name='ff_comb',
        formula='{frac_tt}*{ff_tt} + ({frac_w}+{frac_dy})*{ff_w} + {frac_qcd}*{ff_qcd_os}',
            leaves=[
                # Individual fake factors
                qcd_os,
                w,
                tt,
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
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_QCD_MVis_Iso_SS2OS_up',
                         vars=['mvis', 'mu_iso']
                     ),
                     comb.find('ff_qcd_os')
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
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_QCD_MVis_Iso_SS2OS_down',
                         vars=['mvis', 'mu_iso']
                     ),
                     comb.find('ff_qcd_os')
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
        comb_w_up = replace_nodes(
            comb,
            {'ff_w':
             Node(
                 name='ff_w_up',
                 formula='(1.+{sys_w_up})*{ff_w}',
                 leaves=[
                     Leaf(
                         name='sys_w_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_W_MVis_MT_up',
                         vars=['mvis', 'mt']
                     ),
                     comb.find('ff_w')
                 ]
             )
            }
        )
        comb_w_down = replace_nodes(
            comb,
            {'ff_w':
             Node(
                 name='ff_w_down',
                 formula='max(0.,1.-{sys_w_down})*{ff_w}',
                 leaves=[
                     Leaf(
                         name='sys_w_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/uncertainties_QCD_W.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_W_MVis_MT_down',
                         vars=['mvis', 'mt']
                     ),
                     comb.find('ff_w')
                 ]
             )
            }
        )
        comb_w_up_stat = replace_nodes(
            comb, 
            {'ff_w':
             Node(
                 name='ff_w_up_stat',
                 formula='(1.+{stat_w_up})*{ff_w}',
                 leaves=[
                     Leaf(
                         name='stat_w_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_W_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_HighMT_3D_Iso_Medium_InvertIso_Medium_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     comb.find('ff_w')
                 ]
             )
            }
        )
        comb_w_down_stat = replace_nodes(
            comb, 
            {'ff_w':
             Node(
                 name='ff_w_down_stat',
                 formula='max(0.,1.-{stat_w_down})*{ff_w}',
                 leaves=[
                     Leaf(
                         name='stat_w_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_W_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_HighMT_3D_Iso_Medium_InvertIso_Medium_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     comb.find('ff_w')
                 ]
             )
            }
        )
        comb_tt_up = replace_nodes(
            comb,
            {'ff_tt':
             Node(
                 name='ff_tt_up',
                 formula='(1.+{sys_tt_up})*{ff_tt}',
                 leaves=[
                     Leaf(
                         name='sys_tt_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/uncertainties_TT.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_TT_MVis_up',
                         vars=['tau_decay','mvis']
                     ),
                     comb.find('ff_tt')
                 ]
             )
            }
        )
        comb_tt_down = replace_nodes(
            comb,
            {'ff_tt':
             Node(
                 name='ff_tt_down',
                 formula='max(0.,1.-{sys_tt_down})*{ff_tt}',
                 leaves=[
                     Leaf(
                         name='sys_tt_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/uncertainties_TT.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='uncertainties_TT_MVis_down',
                         vars=['tau_decay','mvis']
                     ),
                     comb.find('ff_tt')
                 ]
             )
            }
        )
        comb_tt_up_stat = replace_nodes(
            comb, 
            {'ff_tt':
             Node(
                 name='ff_tt_up_stat',
                 formula='(1.+{stat_tt_up})*{ff_tt}',
                 leaves=[
                     Leaf(
                         name='stat_tt_up',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_TT_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_TT_anyb_addLep_InvertIso_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     tt.find('ff_tt')
                 ]
             )
            }
        )
        comb_tt_down_stat = replace_nodes(
            comb, 
            {'ff_tt':
             Node(
                 name='ff_tt_down_stat',
                 formula='max(0.,1.-{stat_tt_down})*{ff_tt}',
                 leaves=[
                     Leaf(
                         name='stat_tt_down',
                         file='{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/pieces/FakeFactors_Data_TT_3D.root'.format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category),
                         object='FakeFactors_Data_TT_anyb_addLep_InvertIso_tau_pt_vs_decayMode_error',
                         vars=['tau_pt','tau_decay','njets']
                     ),
                     tt.find('ff_tt')
                 ]
             )
            }
        )
        
        
        
        
        fill(ff_qcd_os, qcd_os)
        fill(ff_qcd_os, qcd_os_up,   sys='ff_qcd_syst_up')
        fill(ff_qcd_os, qcd_os_down, sys='ff_qcd_syst_down')
        fill(ff_qcd_os, qcd_os_up_stat,   sys='ff_qcd_stat_up')
        fill(ff_qcd_os, qcd_os_down_stat, sys='ff_qcd_stat_down')
        fill(ff_w     , w)
        fill(ff_w, w_up,   sys='ff_w_syst_up')
        fill(ff_w, w_down, sys='ff_w_syst_down')
        fill(ff_w, w_up_stat,   sys='ff_w_stat_up')
        fill(ff_w, w_down_stat, sys='ff_w_stat_down')
        fill(ff_tt    , tt)
        fill(ff_tt, tt_up,   sys='ff_tt_syst_up')
        fill(ff_tt, tt_down, sys='ff_tt_syst_down')
        fill(ff_tt, tt_up_stat,   sys='ff_tt_stat_up')
        fill(ff_tt, tt_down_stat, sys='ff_tt_stat_down')
        fill(ff_comb  , comb)
        fill(ff_comb, comb_qcd_up,   sys='ff_qcd_syst_up')
        fill(ff_comb, comb_qcd_down,   sys='ff_qcd_syst_down')
        fill(ff_comb, comb_qcd_up_stat,   sys='ff_qcd_stat_up')
        fill(ff_comb, comb_qcd_down_stat,   sys='ff_qcd_stat_down')
        fill(ff_comb, comb_w_up,   sys='ff_w_syst_up')
        fill(ff_comb, comb_w_down,   sys='ff_w_syst_down')
        fill(ff_comb, comb_w_up_stat,   sys='ff_w_stat_up')
        fill(ff_comb, comb_w_down_stat,   sys='ff_w_stat_down')
        fill(ff_comb, comb_tt_up,   sys='ff_tt_syst_up')
        fill(ff_comb, comb_tt_down, sys='ff_tt_syst_down')
        fill(ff_comb, comb_tt_up_stat,   sys='ff_tt_stat_up')
        fill(ff_comb, comb_tt_down_stat, sys='ff_tt_stat_down')
        
        
        
        file = ROOT.TFile.Open("{HOME}/public/Htautau/FakeRate/{VERSION}/{CHANNEL}/{CATEGORY}/fakeFactors_{VERSION}.root".format(HOME=home,VERSION=version,CHANNEL=channel,CATEGORY=category), "recreate")
        # Write meta-data
        version_ts = ROOT.TString(version)
        tag_ts     = ROOT.TString(tag)
        file.WriteObject(version_ts , "version")
        file.WriteObject(tag_ts     , "tag")
        # Write fake factors
        file.WriteObject(ff_qcd_os.fakefactor  , "ff_qcd_os")
        file.WriteObject(ff_w.fakefactor       , "ff_w")
        file.WriteObject(ff_tt.fakefactor      , "ff_tt")
        file.WriteObject(ff_comb.fakefactor    , "ff_comb")
        #
        file.Close()
    
