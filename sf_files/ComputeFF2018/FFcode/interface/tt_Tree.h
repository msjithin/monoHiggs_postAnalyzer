  float genweight,genM,genpT,pt_top1,pt_top2;

float byTightDPF_2,byVVVLooseDeepVSjet_2,byVVLooseDeepVSjet_2,byVLooseDeepVSjet_2,byLooseDeepVSjet_2,byMediumDeepVSjet_2,byTightDeepVSjet_2,byVTightDeepVSjet_2,byVVTightDeepVSjet_2;
float byVVVLooseDeepVSmu_2,byVVLooseDeepVSmu_2,byVLooseDeepVSmu_2,byLooseDeepVSmu_2,byMediumDeepVSmu_2,byTightDeepVSmu_2,byVTightDeepVSmu_2,byVVTightDeepVSmu_2;
float byVVVLooseDeepVSe_2,byVVLooseDeepVSe_2,byVLooseDeepVSe_2,byLooseDeepVSe_2,byMediumDeepVSe_2,byTightDeepVSe_2,byVTightDeepVSe_2,byVVTightDeepVSe_2;

float byTightDPF_1,byVVVLooseDeepVSjet_1,byVVLooseDeepVSjet_1,byVLooseDeepVSjet_1,byLooseDeepVSjet_1,byMediumDeepVSjet_1,byTightDeepVSjet_1,byVTightDeepVSjet_1,byVVTightDeepVSjet_1;
float byVVVLooseDeepVSmu_1,byVVLooseDeepVSmu_1,byVLooseDeepVSmu_1,byLooseDeepVSmu_1,byMediumDeepVSmu_1,byTightDeepVSmu_1,byVTightDeepVSmu_1,byVVTightDeepVSmu_1;
float byVVVLooseDeepVSe_1,byVVLooseDeepVSe_1,byVLooseDeepVSe_1,byLooseDeepVSe_1,byMediumDeepVSe_1,byTightDeepVSe_1,byVTightDeepVSe_1,byVVTightDeepVSe_1;

float passDoubleTightTau35TightID, passDoubleMediumTau40TightID, passDoubleTightTau40ID, passDoubleMediumHPSTau35, passDoubleTau35, passDoubleTauCmb35;

  float matchEmbFilter_Ele24Tau30_1,matchEmbFilter_Ele27_1,matchEmbFilter_Ele32DoubleL1v1_1,matchEmbFilter_Ele32DoubleL1v2_1,matchEmbFilter_Ele32_1,matchEmbFilter_Ele35_1,matchEmbFilter_Ele24Tau30_2;
  float passEle27, passEle32, passEle35, passEle24Tau30, passEle24HPSTau30, passEle25, matchEle25_1, filterEle25_1;
  float matchEle27_1, matchEle32_1, matchEle35_1, matchEle24Tau30_1, matchEle24Tau30_2, matchEle24HPSTau30_1, matchEle24HPSTau30_2;
  float filterEle27_1, filterEle32_1, filterEle35_1, filterEle24Tau30_1, filterEle24Tau30_2, filterEle24HPSTau30_1, filterEle24HPSTau30_2;
  float photonIso_2, puIso_2, chargedIso_2, neutralIso_2,metcor,metcorphi;
  float byCombinedIsolationDeltaBetaCorrRaw3Hits_2, byIsolationMVA3oldDMwLTraw_2;
  float pt_1,pt_2,px_1,px_2,py_1,py_2,pz_1,pz_2,eta_1,eta_2,phi_1,phi_2,iso_1,e_1,e_2,m_1,m_2,isoDB_1; 
  float pt_1_ScaleUp, pt_1_ScaleDown, pt_1_SigmaUp, pt_1_SigmaDown;
  float flag_BadChargedCandidate, flag_BadPFMuon, flag_EcalDeadCellTriggerPrimitive, flag_HBHENoise, flag_HBHENoiseIso, flag_badCloneMuon, flag_badGlobalMuon, flag_eeBadSc, flag_globalTightHalo2016, flag_goodVertices, flag_globalSuperTightHalo2016, flag_badMuons, flag_duplicateMuons, flag_ecalBadCalib;
  float Flag_BadChargedCandidateFilter, Flag_BadPFMuonFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_HBHENoiseFilter, Flag_HBHENoiseIsoFilter, Flag_badCloneMuon, Flag_badGlobalMuon, Flag_eeBadScFilter, Flag_goodVertices, Flag_globalSuperTightHalo2016Filter, Flag_badMuons, Flag_duplicateMuons, Flag_ecalBadCalibFilter, Flag_ecalBadCalibReducedMINIAODFilter;
  int gen_match_1, gen_match_2;
  int nbtag, nbtagL;
  float bweight;
  float numGenJets;
  int run, lumi, evt;
  float met_JetEta0to3Up, met_JetEta0to5Up, met_JetRelativeSampleUp, met_JetRelativeBalUp, met_JetEta3to5Up, met_JetEta0to3Down, met_JetEta0to5Down, met_JetRelativeSampleDown, met_JetRelativeBalDown, met_JetEta3to5Down, met_JetEC2Up, met_JetEC2Down;
  float metphi_JetEta0to3Up, metphi_JetEta0to5Up, metphi_JetRelativeSampleUp, metphi_JetRelativeBalUp, metphi_JetEta3to5Up, metphi_JetEta0to3Down, metphi_JetEta0to5Down, metphi_JetRelativeSampleDown, metphi_JetRelativeBalDown, metphi_JetEta3to5Down, metphi_JetEC2Up, metphi_JetEC2Down;
  float met_UESUp, met_UESDown, metphi_UESUp, metphi_UESDown, met_responseUp, metphi_responseUp, met_responseDown, metphi_responseDown, met_resolutionUp, metphi_resolutionUp, met_resolutionDown, metphi_resolutionDown;
  float Rivet_stage1_cat_pTjet30GeV,Rivet_stage0_cat, Rivet_higgsPt, Rivet_nJets30, Rivet_stage1p1_cat;
  int nup,njets;
  int njets_JetEta0to3Up, njets_JetEta0to5Up, njets_JetRelativeSampleUp, njets_JetRelativeBalUp, njets_JetEta3to5Up, njets_JetEta0to3Down, njets_JetEta0to5Down, njets_JetRelativeSampleDown, njets_JetRelativeBalDown, njets_JetEta3to5Down, njets_JetEC2Up, njets_JetEC2Down;
  float mjj_JetEta0to3Up, mjj_JetEta0to5Up, mjj_JetRelativeSampleUp, mjj_JetRelativeBalUp, mjj_JetEta3to5Up, mjj_JetEta0to3Down, mjj_JetEta0to5Down, mjj_JetRelativeSampleDown, mjj_JetRelativeBalDown, mjj_JetEta3to5Down, mjj_JetEC2Up, mjj_JetEC2Down;
  float m_sv, m_sv_UP, m_sv_DOWN, m_sv_UESUp, m_sv_UESDown, m_sv_ResolutionUp, m_sv_ResolutionDown, m_sv_ResponseUp, m_sv_ResponseDown, m_sv_JetRelativeSampleUp, m_sv_JetRelativeSampleDown, m_sv_JetRelativeBalUp, m_sv_JetRelativeBalDown, m_sv_JetEta0to3Up, m_sv_JetEta0to3Down, m_sv_JetEta0to5Up, m_sv_JetEta0to5Down, m_sv_JetEta3to5Up, m_sv_JetEta3to5Down, m_sv_JetEC2Up, m_sv_JetEC2Down,m_sv_ESCALEUP, m_sv_ESCALEDOWN, m_sv_ESMEARUP, m_sv_ESMEARDOWN;
  float npv, jpt_1,mjj, jeta_1, jeta_2, jpt_2, jphi_1, jphi_2;
  float q_1,q_2;
  float met, metphi, met_px, met_py;
  float bpt_1, beta_1, bphi_1, bflavor_1, bpt_2, beta_2, bphi_2, bflavor_2;
  float againstElectronVLooseMVA6_2,againstElectronLooseMVA6_2,againstElectronMediumMVA6_2,againstElectronTightMVA6_2,againstElectronVTightMVA6_2;
  float againstElectronVLooseMVA62018_2,againstElectronLooseMVA62018_2,againstElectronMediumMVA62018_2,againstElectronTightMVA62018_2,againstElectronVTightMVA62018_2;
  float againstMuonLoose3_2,againstMuonTight3_2;
  float byVVLooseIsolationMVArun2v2DBoldDMwLT_2, byVLooseIsolationMVArun2v2DBoldDMwLT_2, byLooseIsolationMVArun2v2DBoldDMwLT_2, byMediumIsolationMVArun2v2DBoldDMwLT_2, byTightIsolationMVArun2v2DBoldDMwLT_2, byVTightIsolationMVArun2v2DBoldDMwLT_2, byVVTightIsolationMVArun2v2DBoldDMwLT_2;
  float l2_decayMode, l1_decayMode;
  float decayModeFinding_2;
  float puweight;
  float npu;
  int NUP;

