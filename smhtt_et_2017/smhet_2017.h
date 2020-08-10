//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 23 10:52:42 2020 by ROOT version 6.12/07
// from TTree etau_tree/etau_tree
// found on file: DY_v1.root
//////////////////////////////////////////////////////////

#ifndef smhet_2017_h
#define smhet_2017_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <list>
#include <vector>
#include <bitset>
#include <TCanvas.h>
#include <TSystem.h>
#include <TPostScript.h>
#include <TH2.h>
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TRef.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include"RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"
#include "vector"
#include "TString.h"
#include "vector"

using namespace std;
// Header file for the classes stored in the TTree if any.

class smhet_2017 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   TFile *fileName;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Float_t         genpX;
   Float_t         genpY;
   Float_t         genM;
   Float_t         genpT;
   Float_t         vispX;
   Float_t         vispY;
   Float_t         Rivet_VEta;
   Float_t         Rivet_VPt;
   Float_t         Rivet_errorCode;
   Float_t         Rivet_higgsEta;
   Float_t         Rivet_higgsPt;
   Float_t         Rivet_nJets25;
   Float_t         Rivet_nJets30;
   Float_t         Rivet_p4decay_VEta;
   Float_t         Rivet_p4decay_VPt;
   Float_t         Rivet_prodMode;
   Float_t         Rivet_stage0_cat;
   Float_t         Rivet_stage1p1_cat;
   Float_t         Rivet_stage1_cat_pTjet25GeV;
   Float_t         Rivet_stage1_cat_pTjet30GeV;
   Float_t         npv;
   Float_t         npu;
   Float_t         L1iso;
   Float_t         L1pt;
   Float_t         pt_1;
   Float_t         pt_1_ScaleUp;
   Float_t         pt_1_ScaleDown;
   Float_t         pt_1_SigmaUp;
   Float_t         pt_1_SigmaDown;
   Float_t         phi_1;
   Float_t         eta_1;
   Float_t         m_1;
   Float_t         e_1;
   Float_t         q_1;
   Float_t         iso_1;
   Float_t         isoDB_1;
   Float_t         eid80_iso_1;
   Float_t         eid90_iso_1;
   Float_t         eid80_noiso_1;
   Float_t         eid90_noiso_1;
   Float_t         genpt_1;
   Float_t         geneta_1;
   Float_t         genpt_2;
   Float_t         geneta_2;
   Float_t         prefiring_weight;
   Float_t         prefiring_weight_up;
   Float_t         prefiring_weight_down;
   Float_t         bweight;
   Float_t         pt_2;
   Float_t         phi_2;
   Float_t         eta_2;
   Float_t         m_2;
   Float_t         e_2;
   Float_t         q_2;
   Float_t         l2_decayMode;
   Float_t         againstElectronVLooseMVA6_2;
   Float_t         againstElectronLooseMVA6_2;
   Float_t         againstElectronMediumMVA6_2;
   Float_t         againstElectronTightMVA6_2;
   Float_t         againstElectronVTightMVA6_2;
   Float_t         againstElectronVLooseMVA62018_2;
   Float_t         againstElectronLooseMVA62018_2;
   Float_t         againstElectronMediumMVA62018_2;
   Float_t         againstElectronTightMVA62018_2;
   Float_t         againstElectronVTightMVA62018_2;
   Float_t         againstMuonLoose3_2;
   Float_t         againstMuonTight3_2;
   Float_t         byVVLooseIsolationMVArun2v2DBoldDMwLT_2;
   Float_t         byVLooseIsolationMVArun2v2DBoldDMwLT_2;
   Float_t         byLooseIsolationMVArun2v2DBoldDMwLT_2;
   Float_t         byMediumIsolationMVArun2v2DBoldDMwLT_2;
   Float_t         byTightIsolationMVArun2v2DBoldDMwLT_2;
   Float_t         byVTightIsolationMVArun2v2DBoldDMwLT_2;
   Float_t         byVVTightIsolationMVArun2v2DBoldDMwLT_2;
   Float_t         byIsolationMVA3oldDMwLTraw_2;
   Float_t         decayModeFinding_2;
   Float_t         byTightDPF_2;
   Float_t         byVVVLooseDeepVSjet_2;
   Float_t         byVVLooseDeepVSjet_2;
   Float_t         byVLooseDeepVSjet_2;
   Float_t         byLooseDeepVSjet_2;
   Float_t         byMediumDeepVSjet_2;
   Float_t         byTightDeepVSjet_2;
   Float_t         byVTightDeepVSjet_2;
   Float_t         byVVTightDeepVSjet_2;
   Float_t         byVVVLooseDeepVSmu_2;
   Float_t         byVVLooseDeepVSmu_2;
   Float_t         byVLooseDeepVSmu_2;
   Float_t         byLooseDeepVSmu_2;
   Float_t         byMediumDeepVSmu_2;
   Float_t         byTightDeepVSmu_2;
   Float_t         byVTightDeepVSmu_2;
   Float_t         byVVTightDeepVSmu_2;
   Float_t         byVVVLooseDeepVSe_2;
   Float_t         byVVLooseDeepVSe_2;
   Float_t         byVLooseDeepVSe_2;
   Float_t         byLooseDeepVSe_2;
   Float_t         byMediumDeepVSe_2;
   Float_t         byTightDeepVSe_2;
   Float_t         byVTightDeepVSe_2;
   Float_t         byVVTightDeepVSe_2;
   Float_t         numGenJets;
   Float_t         jetPt_2;
   Float_t         Flag_goodVertices;
   Float_t         Flag_globalSuperTightHalo2016Filter;
   Float_t         Flag_eeBadScFilter;
   Float_t         Flag_ecalBadCalibFilter;
   Float_t         Flag_ecalBadCalibReducedMINIAODFilter;
   Float_t         Flag_badMuons;
   Float_t         Flag_duplicateMuons;
   Float_t         Flag_HBHENoiseIsoFilter;
   Float_t         Flag_HBHENoiseFilter;
   Float_t         Flag_EcalDeadCellTriggerPrimitiveFilter;
   Float_t         Flag_BadPFMuonFilter;
   Float_t         Flag_BadChargedCandidateFilter;
   Float_t         matchEmbFilter_Ele24Tau30_2;
   Float_t         matchEmbFilter_Ele27_1;
   Float_t         matchEmbFilter_Ele32DoubleL1v2_1;
   Float_t         matchEmbFilter_Ele32DoubleL1v1_1;
   Float_t         matchEmbFilter_Ele32_1;
   Float_t         matchEmbFilter_Ele35_1;
   Float_t         matchEmbFilter_Ele24Tau30_1;
   Float_t         met;
   Float_t         metSig;
   Float_t         metcov00;
   Float_t         metcov10;
   Float_t         metcov11;
   Float_t         metcov01;
   Float_t         metphi;
   Float_t         met_py;
   Float_t         met_px;
   Float_t         met_UESUp;
   Float_t         metphi_UESUp;
   Float_t         met_UESDown;
   Float_t         metphi_UESDown;
   Float_t         met_responseUp;
   Float_t         met_responseDown;
   Float_t         met_resolutionUp;
   Float_t         met_resolutionDown;
   Float_t         metphi_responseUp;
   Float_t         metphi_responseDown;
   Float_t         metphi_resolutionUp;
   Float_t         metphi_resolutionDown;
   Float_t         met_JetEta0to3Up;
   Float_t         metphi_JetEta0to3Up;
   Float_t         met_JetEta0to3Down;
   Float_t         metphi_JetEta0to3Down;
   Float_t         met_JetEC2Up;
   Float_t         metphi_JetEC2Up;
   Float_t         met_JetEC2Down;
   Float_t         metphi_JetEC2Down;
   Float_t         met_JetRelativeBalUp;
   Float_t         metphi_JetRelativeBalUp;
   Float_t         met_JetRelativeBalDown;
   Float_t         metphi_JetRelativeBalDown;
   Float_t         met_JetRelativeSampleUp;
   Float_t         metphi_JetRelativeSampleUp;
   Float_t         met_JetRelativeSampleDown;
   Float_t         metphi_JetRelativeSampleDown;
   Float_t         met_JetEta3to5Up;
   Float_t         metphi_JetEta3to5Up;
   Float_t         met_JetEta3to5Down;
   Float_t         metphi_JetEta3to5Down;
   Float_t         met_JetEta0to5Up;
   Float_t         metphi_JetEta0to5Up;
   Float_t         met_JetEta0to5Down;
   Float_t         metphi_JetEta0to5Down;
   Float_t         passEle27;
   Float_t         passEle32;
   Float_t         passEle35;
   Float_t         passEle24Tau30;
   Float_t         matchEle27_1;
   Float_t         matchEle32_1;
   Float_t         matchEle35_1;
   Float_t         matchEle24Tau30_1;
   Float_t         matchEle24Tau30_2;
   Float_t         filterEle27_1;
   Float_t         filterEle32_1;
   Float_t         filterEle35_1;
   Float_t         filterEle24Tau30_1;
   Float_t         filterEle24Tau30_2;
   Float_t         mjj;
   Float_t         mjj_JetRelativeBalUp;
   Float_t         mjj_JetRelativeBalDown;
   Float_t         mjj_JetRelativeSampleUp;
   Float_t         mjj_JetRelativeSampleDown;
   Float_t         mjj_JetEta0to3Up;
   Float_t         mjj_JetEta0to3Down;
   Float_t         mjj_JetEC2Up;
   Float_t         mjj_JetEC2Down;
   Float_t         mjj_JetEta0to5Up;
   Float_t         mjj_JetEta0to5Down;
   Float_t         mjj_JetEta3to5Up;
   Float_t         mjj_JetEta3to5Down;
   Int_t           gen_match_1;
   Int_t           gen_match_2;
   Int_t           nbtag;
   Int_t           nbtagL;
   Int_t           njets;
   Int_t           njets_JetRelativeBalUp;
   Int_t           njets_JetRelativeBalDown;
   Int_t           njets_JetRelativeSampleUp;
   Int_t           njets_JetRelativeSampleDown;
   Int_t           njets_JetEta0to3Up;
   Int_t           njets_JetEta0to3Down;
   Int_t           njets_JetEC2Up;
   Int_t           njets_JetEC2Down;
   Int_t           njets_JetEta0to5Up;
   Int_t           njets_JetEta0to5Down;
   Int_t           njets_JetEta3to5Up;
   Int_t           njets_JetEta3to5Down;
   Float_t         jpt_1;
   Float_t         jeta_1;
   Float_t         jphi_1;
   Float_t         jpt_2;
   Float_t         jeta_2;
   Float_t         jphi_2;
   Float_t         jcsv_1;
   Float_t         jcsv_2;
   Float_t         jpt_JetEta0to3Up_1;
   Float_t         jpt_JetEta0to3Down_1;
   Float_t         jpt_JetEC2Up_1;
   Float_t         jpt_JetEC2Down_1;
   Float_t         jpt_JetEta3to5Up_1;
   Float_t         jpt_JetEta3to5Down_1;
   Float_t         jpt_JetEta0to5Up_1;
   Float_t         jpt_JetEta0to5Down_1;
   Float_t         jpt_JetRelativeBalUp_1;
   Float_t         jpt_JetRelativeBalDown_1;
   Float_t         jpt_JetRelativeSampleUp_1;
   Float_t         jpt_JetRelativeSampleDown_1;
   Float_t         jpt_JetEta0to3Up_2;
   Float_t         jpt_JetEta0to3Down_2;
   Float_t         jpt_JetEC2Up_2;
   Float_t         jpt_JetEC2Down_2;
   Float_t         jpt_JetEta3to5Up_2;
   Float_t         jpt_JetEta3to5Down_2;
   Float_t         jpt_JetEta0to5Up_2;
   Float_t         jpt_JetEta0to5Down_2;
   Float_t         jpt_JetRelativeBalUp_2;
   Float_t         jpt_JetRelativeBalDown_2;
   Float_t         jpt_JetRelativeSampleUp_2;
   Float_t         jpt_JetRelativeSampleDown_2;
   Float_t         bpt_1;
   Float_t         bflavor_1;
   Float_t         beta_1;
   Float_t         bphi_1;
   Float_t         bpt_2;
   Float_t         bflavor_2;
   Float_t         beta_2;
   Float_t         bphi_2;
   Float_t         pt_top1;
   Float_t         pt_top2;
   Float_t         genweight;
   Float_t         gen_Higgs_pt;
   Float_t         gen_Higgs_mass;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_genpX;   //!
   TBranch        *b_genpY;   //!
   TBranch        *b_genM;   //!
   TBranch        *b_genpT;   //!
   TBranch        *b_vispX;   //!
   TBranch        *b_vispY;   //!
   TBranch        *b_Rivet_VEta;   //!
   TBranch        *b_Rivet_VPt;   //!
   TBranch        *b_Rivet_errorCode;   //!
   TBranch        *b_Rivet_higgsEta;   //!
   TBranch        *b_Rivet_higgsPt;   //!
   TBranch        *b_Rivet_nJets25;   //!
   TBranch        *b_Rivet_nJets30;   //!
   TBranch        *b_Rivet_p4decay_VEta;   //!
   TBranch        *b_Rivet_p4decay_VPt;   //!
   TBranch        *b_Rivet_prodMode;   //!
   TBranch        *b_Rivet_stage0_cat;   //!
   TBranch        *b_Rivet_stage1p1_cat;   //!
   TBranch        *b_Rivet_stage1_cat_pTjet25GeV;   //!
   TBranch        *b_Rivet_stage1_cat_pTjet30GeV;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_L1iso;   //!
   TBranch        *b_L1pt;   //!
   TBranch        *b_pt_1;   //!
   TBranch        *b_pt_1_ScaleUp;   //!
   TBranch        *b_pt_1_ScaleDown;   //!
   TBranch        *b_pt_1_SigmaUp;   //!
   TBranch        *b_pt_1_SigmaDown;   //!
   TBranch        *b_phi_1;   //!
   TBranch        *b_eta_1;   //!
   TBranch        *b_m_1;   //!
   TBranch        *b_e_1;   //!
   TBranch        *b_q_1;   //!
   TBranch        *b_iso_1;   //!
   TBranch        *b_isoDB_1;   //!
   TBranch        *b_eid80_iso_1;   //!
   TBranch        *b_eid90_iso_1;   //!
   TBranch        *b_eid80_noiso_1;   //!
   TBranch        *b_eid90_noiso_1;   //!
   TBranch        *b_genpt_1;   //!
   TBranch        *b_geneta_1;   //!
   TBranch        *b_genpt_2;   //!
   TBranch        *b_geneta_2;   //!
   TBranch        *b_prefiring_weight;   //!
   TBranch        *b_prefiring_weight_up;   //!
   TBranch        *b_prefiring_weight_down;   //!
   TBranch        *b_bweight;   //!
   TBranch        *b_pt_2;   //!
   TBranch        *b_phi_2;   //!
   TBranch        *b_eta_2;   //!
   TBranch        *b_m_2;   //!
   TBranch        *b_e_2;   //!
   TBranch        *b_q_2;   //!
   TBranch        *b_l2_decayMode;   //!
   TBranch        *b_againstElectronVLooseMVA6_2;   //!
   TBranch        *b_againstElectronLooseMVA6_2;   //!
   TBranch        *b_againstElectronMediumMVA6_2;   //!
   TBranch        *b_againstElectronTightMVA6_2;   //!
   TBranch        *b_againstElectronVTightMVA6_2;   //!
   TBranch        *b_againstElectronVLooseMVA62018_2;   //!
   TBranch        *b_againstElectronLooseMVA62018_2;   //!
   TBranch        *b_againstElectronMediumMVA62018_2;   //!
   TBranch        *b_againstElectronTightMVA62018_2;   //!
   TBranch        *b_againstElectronVTightMVA62018_2;   //!
   TBranch        *b_againstMuonLoose3_2;   //!
   TBranch        *b_againstMuonTight3_2;   //!
   TBranch        *b_byVVLooseIsolationMVArun2v2DBoldDMwLT_2;   //!
   TBranch        *b_byVLooseIsolationMVArun2v2DBoldDMwLT_2;   //!
   TBranch        *b_byLooseIsolationMVArun2v2DBoldDMwLT_2;   //!
   TBranch        *b_byMediumIsolationMVArun2v2DBoldDMwLT_2;   //!
   TBranch        *b_byTightIsolationMVArun2v2DBoldDMwLT_2;   //!
   TBranch        *b_byVTightIsolationMVArun2v2DBoldDMwLT_2;   //!
   TBranch        *b_byVVTightIsolationMVArun2v2DBoldDMwLT_2;   //!
   TBranch        *b_byIsolationMVA3oldDMwLTraw_2;   //!
   TBranch        *b_decayModeFinding_2;   //!
   TBranch        *b_byTightDPF_2;   //!
   TBranch        *b_byVVVLooseDeepVSjet_2;   //!
   TBranch        *b_byVVLooseDeepVSjet_2;   //!
   TBranch        *b_byVLooseDeepVSjet_2;   //!
   TBranch        *b_byLooseDeepVSjet_2;   //!
   TBranch        *b_byMediumDeepVSjet_2;   //!
   TBranch        *b_byTightDeepVSjet_2;   //!
   TBranch        *b_byVTightDeepVSjet_2;   //!
   TBranch        *b_byVVTightDeepVSjet_2;   //!
   TBranch        *b_byVVVLooseDeepVSmu_2;   //!
   TBranch        *b_byVVLooseDeepVSmu_2;   //!
   TBranch        *b_byVLooseDeepVSmu_2;   //!
   TBranch        *b_byLooseDeepVSmu_2;   //!
   TBranch        *b_byMediumDeepVSmu_2;   //!
   TBranch        *b_byTightDeepVSmu_2;   //!
   TBranch        *b_byVTightDeepVSmu_2;   //!
   TBranch        *b_byVVTightDeepVSmu_2;   //!
   TBranch        *b_byVVVLooseDeepVSe_2;   //!
   TBranch        *b_byVVLooseDeepVSe_2;   //!
   TBranch        *b_byVLooseDeepVSe_2;   //!
   TBranch        *b_byLooseDeepVSe_2;   //!
   TBranch        *b_byMediumDeepVSe_2;   //!
   TBranch        *b_byTightDeepVSe_2;   //!
   TBranch        *b_byVTightDeepVSe_2;   //!
   TBranch        *b_byVVTightDeepVSe_2;   //!
   TBranch        *b_numGenJets;   //!
   TBranch        *b_jetPt_2;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_ecalBadCalibReducedMINIAODFilter;   //!
   TBranch        *b_Flag_badMuons;   //!
   TBranch        *b_Flag_duplicateMuons;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_matchEmbFilter_Ele24Tau30_2;   //!
   TBranch        *b_matchEmbFilter_Ele27_1;   //!
   TBranch        *b_matchEmbFilter_Ele32DoubleL1v2_1;   //!
   TBranch        *b_matchEmbFilter_Ele32DoubleL1v1_1;   //!
   TBranch        *b_matchEmbFilter_Ele32_1;   //!
   TBranch        *b_matchEmbFilter_Ele35_1;   //!
   TBranch        *b_matchEmbFilter_Ele24Tau30_1;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metSig;   //!
   TBranch        *b_metcov00;   //!
   TBranch        *b_metcov10;   //!
   TBranch        *b_metcov11;   //!
   TBranch        *b_metcov01;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_met_py;   //!
   TBranch        *b_met_px;   //!
   TBranch        *b_met_UESUp;   //!
   TBranch        *b_metphi_UESUp;   //!
   TBranch        *b_met_UESDown;   //!
   TBranch        *b_metphi_UESDown;   //!
   TBranch        *b_met_responseUp;   //!
   TBranch        *b_met_responseDown;   //!
   TBranch        *b_met_resolutionUp;   //!
   TBranch        *b_met_resolutionDown;   //!
   TBranch        *b_metphi_responseUp;   //!
   TBranch        *b_metphi_responseDown;   //!
   TBranch        *b_metphi_resolutionUp;   //!
   TBranch        *b_metphi_resolutionDown;   //!
   TBranch        *b_met_JetEta0to3Up;   //!
   TBranch        *b_metphi_JetEta0to3Up;   //!
   TBranch        *b_met_JetEta0to3Down;   //!
   TBranch        *b_metphi_JetEta0to3Down;   //!
   TBranch        *b_met_JetEC2Up;   //!
   TBranch        *b_metphi_JetEC2Up;   //!
   TBranch        *b_met_JetEC2Down;   //!
   TBranch        *b_metphi_JetEC2Down;   //!
   TBranch        *b_met_JetRelativeBalUp;   //!
   TBranch        *b_metphi_JetRelativeBalUp;   //!
   TBranch        *b_met_JetRelativeBalDown;   //!
   TBranch        *b_metphi_JetRelativeBalDown;   //!
   TBranch        *b_met_JetRelativeSampleUp;   //!
   TBranch        *b_metphi_JetRelativeSampleUp;   //!
   TBranch        *b_met_JetRelativeSampleDown;   //!
   TBranch        *b_metphi_JetRelativeSampleDown;   //!
   TBranch        *b_met_JetEta3to5Up;   //!
   TBranch        *b_metphi_JetEta3to5Up;   //!
   TBranch        *b_met_JetEta3to5Down;   //!
   TBranch        *b_metphi_JetEta3to5Down;   //!
   TBranch        *b_met_JetEta0to5Up;   //!
   TBranch        *b_metphi_JetEta0to5Up;   //!
   TBranch        *b_met_JetEta0to5Down;   //!
   TBranch        *b_metphi_JetEta0to5Down;   //!
   TBranch        *b_passEle27;   //!
   TBranch        *b_passEle32;   //!
   TBranch        *b_passEle35;   //!
   TBranch        *b_passEle24Tau30;   //!
   TBranch        *b_matchEle27_1;   //!
   TBranch        *b_matchEle32_1;   //!
   TBranch        *b_matchEle35_1;   //!
   TBranch        *b_matchEle24Tau30_1;   //!
   TBranch        *b_matchEle24Tau30_2;   //!
   TBranch        *b_filterEle27_1;   //!
   TBranch        *b_filterEle32_1;   //!
   TBranch        *b_filterEle35_1;   //!
   TBranch        *b_filterEle24Tau30_1;   //!
   TBranch        *b_filterEle24Tau30_2;   //!
   TBranch        *b_mjj;   //!
   TBranch        *b_mjj_JetRelativeBalUp;   //!
   TBranch        *b_mjj_JetRelativeBalDown;   //!
   TBranch        *b_mjj_JetRelativeSampleUp;   //!
   TBranch        *b_mjj_JetRelativeSampleDown;   //!
   TBranch        *b_mjj_JetEta0to3Up;   //!
   TBranch        *b_mjj_JetEta0to3Down;   //!
   TBranch        *b_mjj_JetEC2Up;   //!
   TBranch        *b_mjj_JetEC2Down;   //!
   TBranch        *b_mjj_JetEta0to5Up;   //!
   TBranch        *b_mjj_JetEta0to5Down;   //!
   TBranch        *b_mjj_JetEta3to5Up;   //!
   TBranch        *b_mjj_JetEta3to5Down;   //!
   TBranch        *b_gen_match_1;   //!
   TBranch        *b_gen_match_2;   //!
   TBranch        *b_nbtag;   //!
   TBranch        *b_nbtagL;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_njets_JetRelativeBalUp;   //!
   TBranch        *b_njets_JetRelativeBalDown;   //!
   TBranch        *b_njets_JetRelativeSampleUp;   //!
   TBranch        *b_njets_JetRelativeSampleDown;   //!
   TBranch        *b_njets_JetEta0to3Up;   //!
   TBranch        *b_njets_JetEta0to3Down;   //!
   TBranch        *b_njets_JetEC2Up;   //!
   TBranch        *b_njets_JetEC2Down;   //!
   TBranch        *b_njets_JetEta0to5Up;   //!
   TBranch        *b_njets_JetEta0to5Down;   //!
   TBranch        *b_njets_JetEta3to5Up;   //!
   TBranch        *b_njets_JetEta3to5Down;   //!
   TBranch        *b_jpt_1;   //!
   TBranch        *b_jeta_1;   //!
   TBranch        *b_jphi_1;   //!
   TBranch        *b_jpt_2;   //!
   TBranch        *b_jeta_2;   //!
   TBranch        *b_jphi_2;   //!
   TBranch        *b_jcsv_1;   //!
   TBranch        *b_jcsv_2;   //!
   TBranch        *b_jpt_JetEta0to3Up_1;   //!
   TBranch        *b_jpt_JetEta0to3Down_1;   //!
   TBranch        *b_jpt_JetEC2Up_1;   //!
   TBranch        *b_jpt_JetEC2Down_1;   //!
   TBranch        *b_jpt_JetEta3to5Up_1;   //!
   TBranch        *b_jpt_JetEta3to5Down_1;   //!
   TBranch        *b_jpt_JetEta0to5Up_1;   //!
   TBranch        *b_jpt_JetEta0to5Down_1;   //!
   TBranch        *b_jpt_JetRelativeBalUp_1;   //!
   TBranch        *b_jpt_JetRelativeBalDown_1;   //!
   TBranch        *b_jpt_JetRelativeSampleUp_1;   //!
   TBranch        *b_jpt_JetRelativeSampleDown_1;   //!
   TBranch        *b_jpt_JetEta0to3Up_2;   //!
   TBranch        *b_jpt_JetEta0to3Down_2;   //!
   TBranch        *b_jpt_JetEC2Up_2;   //!
   TBranch        *b_jpt_JetEC2Down_2;   //!
   TBranch        *b_jpt_JetEta3to5Up_2;   //!
   TBranch        *b_jpt_JetEta3to5Down_2;   //!
   TBranch        *b_jpt_JetEta0to5Up_2;   //!
   TBranch        *b_jpt_JetEta0to5Down_2;   //!
   TBranch        *b_jpt_JetRelativeBalUp_2;   //!
   TBranch        *b_jpt_JetRelativeBalDown_2;   //!
   TBranch        *b_jpt_JetRelativeSampleUp_2;   //!
   TBranch        *b_jpt_JetRelativeSampleDown_2;   //!
   TBranch        *b_bpt_1;   //!
   TBranch        *b_bflavor_1;   //!
   TBranch        *b_beta_1;   //!
   TBranch        *b_bphi_1;   //!
   TBranch        *b_bpt_2;   //!
   TBranch        *b_bflavor_2;   //!
   TBranch        *b_beta_2;   //!
   TBranch        *b_bphi_2;   //!
   TBranch        *b_pt_top1;   //!
   TBranch        *b_pt_top2;   //!
   TBranch        *b_genweight;   //!
   TBranch        *b_gen_Higgs_pt;   //!
   TBranch        *b_gen_Higgs_mass;   //!

   //smhet_2017(TTree *tree=0);
   smhet_2017(const char* file1, const char* file2, string isMC);
   virtual ~smhet_2017();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void BookHistos(const char* file1, const char* file2);
   virtual void     Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_);
   virtual double delta_R(float phi1, float eta1, float phi2, float eta2);
   virtual double DeltaPhi(double phi1, double phi2);
   virtual float TMass_F(float LepPt, float LepPhi , float met, float metPhi);
   virtual void     makeMyPlot( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef smhet_2017_cxx
/* smhet_2017::smhet_2017(TTree *tree) : fChain(0)  */
/* { */
/* // if parameter tree is not specified (or zero), connect the file */
/* // used to generate this class and read the Tree. */
/*    if (tree == 0) { */
/*       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DY_v1.root"); */
/*       if (!f || !f->IsOpen()) { */
/*          f = new TFile("DY_v1.root"); */
/*       } */
/*       f->GetObject("etau_tree",tree); */

/*    } */
/*    Init(tree); */
/* } */
smhet_2017::smhet_2017(const char* file1, const char* file2, string isMC)
{
  TChain *chain = new TChain("etau_tree");
  
  TString  FullPathInputFile = file1;
  chain->Add(FullPathInputFile);
  
  std::cout<<"All files added."<<std::endl;
  std::cout<<"Initializing chain."<<std::endl;
  Init(chain);
  BookHistos(file1, file2);

  //inspected_events->Write();
}

smhet_2017::~smhet_2017()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   fileName->Write();
   fileName->Close();
}

Int_t smhet_2017::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t smhet_2017::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void smhet_2017::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("genpX", &genpX, &b_genpX);
   fChain->SetBranchAddress("genpY", &genpY, &b_genpY);
   fChain->SetBranchAddress("genM", &genM, &b_genM);
   fChain->SetBranchAddress("genpT", &genpT, &b_genpT);
   fChain->SetBranchAddress("vispX", &vispX, &b_vispX);
   fChain->SetBranchAddress("vispY", &vispY, &b_vispY);
   fChain->SetBranchAddress("Rivet_VEta", &Rivet_VEta, &b_Rivet_VEta);
   fChain->SetBranchAddress("Rivet_VPt", &Rivet_VPt, &b_Rivet_VPt);
   fChain->SetBranchAddress("Rivet_errorCode", &Rivet_errorCode, &b_Rivet_errorCode);
   fChain->SetBranchAddress("Rivet_higgsEta", &Rivet_higgsEta, &b_Rivet_higgsEta);
   fChain->SetBranchAddress("Rivet_higgsPt", &Rivet_higgsPt, &b_Rivet_higgsPt);
   fChain->SetBranchAddress("Rivet_nJets25", &Rivet_nJets25, &b_Rivet_nJets25);
   fChain->SetBranchAddress("Rivet_nJets30", &Rivet_nJets30, &b_Rivet_nJets30);
   fChain->SetBranchAddress("Rivet_p4decay_VEta", &Rivet_p4decay_VEta, &b_Rivet_p4decay_VEta);
   fChain->SetBranchAddress("Rivet_p4decay_VPt", &Rivet_p4decay_VPt, &b_Rivet_p4decay_VPt);
   fChain->SetBranchAddress("Rivet_prodMode", &Rivet_prodMode, &b_Rivet_prodMode);
   fChain->SetBranchAddress("Rivet_stage0_cat", &Rivet_stage0_cat, &b_Rivet_stage0_cat);
   fChain->SetBranchAddress("Rivet_stage1p1_cat", &Rivet_stage1p1_cat, &b_Rivet_stage1p1_cat);
   fChain->SetBranchAddress("Rivet_stage1_cat_pTjet25GeV", &Rivet_stage1_cat_pTjet25GeV, &b_Rivet_stage1_cat_pTjet25GeV);
   fChain->SetBranchAddress("Rivet_stage1_cat_pTjet30GeV", &Rivet_stage1_cat_pTjet30GeV, &b_Rivet_stage1_cat_pTjet30GeV);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("L1iso", &L1iso, &b_L1iso);
   fChain->SetBranchAddress("L1pt", &L1pt, &b_L1pt);
   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   fChain->SetBranchAddress("pt_1_ScaleUp", &pt_1_ScaleUp, &b_pt_1_ScaleUp);
   fChain->SetBranchAddress("pt_1_ScaleDown", &pt_1_ScaleDown, &b_pt_1_ScaleDown);
   fChain->SetBranchAddress("pt_1_SigmaUp", &pt_1_SigmaUp, &b_pt_1_SigmaUp);
   fChain->SetBranchAddress("pt_1_SigmaDown", &pt_1_SigmaDown, &b_pt_1_SigmaDown);
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1);
   fChain->SetBranchAddress("m_1", &m_1, &b_m_1);
   fChain->SetBranchAddress("e_1", &e_1, &b_e_1);
   fChain->SetBranchAddress("q_1", &q_1, &b_q_1);
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
   fChain->SetBranchAddress("isoDB_1", &isoDB_1, &b_isoDB_1);
   fChain->SetBranchAddress("eid80_iso_1", &eid80_iso_1, &b_eid80_iso_1);
   fChain->SetBranchAddress("eid90_iso_1", &eid90_iso_1, &b_eid90_iso_1);
   fChain->SetBranchAddress("eid80_noiso_1", &eid80_noiso_1, &b_eid80_noiso_1);
   fChain->SetBranchAddress("eid90_noiso_1", &eid90_noiso_1, &b_eid90_noiso_1);
   fChain->SetBranchAddress("genpt_1", &genpt_1, &b_genpt_1);
   fChain->SetBranchAddress("geneta_1", &geneta_1, &b_geneta_1);
   fChain->SetBranchAddress("genpt_2", &genpt_2, &b_genpt_2);
   fChain->SetBranchAddress("geneta_2", &geneta_2, &b_geneta_2);
   fChain->SetBranchAddress("prefiring_weight", &prefiring_weight, &b_prefiring_weight);
   fChain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up, &b_prefiring_weight_up);
   fChain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down, &b_prefiring_weight_down);
   fChain->SetBranchAddress("bweight", &bweight, &b_bweight);
   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
   fChain->SetBranchAddress("m_2", &m_2, &b_m_2);
   fChain->SetBranchAddress("e_2", &e_2, &b_e_2);
   fChain->SetBranchAddress("q_2", &q_2, &b_q_2);
   fChain->SetBranchAddress("l2_decayMode", &l2_decayMode, &b_l2_decayMode);
   fChain->SetBranchAddress("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2, &b_againstElectronVLooseMVA6_2);
   fChain->SetBranchAddress("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2, &b_againstElectronLooseMVA6_2);
   fChain->SetBranchAddress("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2, &b_againstElectronMediumMVA6_2);
   fChain->SetBranchAddress("againstElectronTightMVA6_2", &againstElectronTightMVA6_2, &b_againstElectronTightMVA6_2);
   fChain->SetBranchAddress("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2, &b_againstElectronVTightMVA6_2);
   fChain->SetBranchAddress("againstElectronVLooseMVA62018_2", &againstElectronVLooseMVA62018_2, &b_againstElectronVLooseMVA62018_2);
   fChain->SetBranchAddress("againstElectronLooseMVA62018_2", &againstElectronLooseMVA62018_2, &b_againstElectronLooseMVA62018_2);
   fChain->SetBranchAddress("againstElectronMediumMVA62018_2", &againstElectronMediumMVA62018_2, &b_againstElectronMediumMVA62018_2);
   fChain->SetBranchAddress("againstElectronTightMVA62018_2", &againstElectronTightMVA62018_2, &b_againstElectronTightMVA62018_2);
   fChain->SetBranchAddress("againstElectronVTightMVA62018_2", &againstElectronVTightMVA62018_2, &b_againstElectronVTightMVA62018_2);
   fChain->SetBranchAddress("againstMuonLoose3_2", &againstMuonLoose3_2, &b_againstMuonLoose3_2);
   fChain->SetBranchAddress("againstMuonTight3_2", &againstMuonTight3_2, &b_againstMuonTight3_2);
   fChain->SetBranchAddress("byVVLooseIsolationMVArun2v2DBoldDMwLT_2", &byVVLooseIsolationMVArun2v2DBoldDMwLT_2, &b_byVVLooseIsolationMVArun2v2DBoldDMwLT_2);
   fChain->SetBranchAddress("byVLooseIsolationMVArun2v2DBoldDMwLT_2", &byVLooseIsolationMVArun2v2DBoldDMwLT_2, &b_byVLooseIsolationMVArun2v2DBoldDMwLT_2);
   fChain->SetBranchAddress("byLooseIsolationMVArun2v2DBoldDMwLT_2", &byLooseIsolationMVArun2v2DBoldDMwLT_2, &b_byLooseIsolationMVArun2v2DBoldDMwLT_2);
   fChain->SetBranchAddress("byMediumIsolationMVArun2v2DBoldDMwLT_2", &byMediumIsolationMVArun2v2DBoldDMwLT_2, &b_byMediumIsolationMVArun2v2DBoldDMwLT_2);
   fChain->SetBranchAddress("byTightIsolationMVArun2v2DBoldDMwLT_2", &byTightIsolationMVArun2v2DBoldDMwLT_2, &b_byTightIsolationMVArun2v2DBoldDMwLT_2);
   fChain->SetBranchAddress("byVTightIsolationMVArun2v2DBoldDMwLT_2", &byVTightIsolationMVArun2v2DBoldDMwLT_2, &b_byVTightIsolationMVArun2v2DBoldDMwLT_2);
   fChain->SetBranchAddress("byVVTightIsolationMVArun2v2DBoldDMwLT_2", &byVVTightIsolationMVArun2v2DBoldDMwLT_2, &b_byVVTightIsolationMVArun2v2DBoldDMwLT_2);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2, &b_byIsolationMVA3oldDMwLTraw_2);
   fChain->SetBranchAddress("decayModeFinding_2", &decayModeFinding_2, &b_decayModeFinding_2);
   fChain->SetBranchAddress("byTightDPF_2", &byTightDPF_2, &b_byTightDPF_2);
   fChain->SetBranchAddress("byVVVLooseDeepVSjet_2", &byVVVLooseDeepVSjet_2, &b_byVVVLooseDeepVSjet_2);
   fChain->SetBranchAddress("byVVLooseDeepVSjet_2", &byVVLooseDeepVSjet_2, &b_byVVLooseDeepVSjet_2);
   fChain->SetBranchAddress("byVLooseDeepVSjet_2", &byVLooseDeepVSjet_2, &b_byVLooseDeepVSjet_2);
   fChain->SetBranchAddress("byLooseDeepVSjet_2", &byLooseDeepVSjet_2, &b_byLooseDeepVSjet_2);
   fChain->SetBranchAddress("byMediumDeepVSjet_2", &byMediumDeepVSjet_2, &b_byMediumDeepVSjet_2);
   fChain->SetBranchAddress("byTightDeepVSjet_2", &byTightDeepVSjet_2, &b_byTightDeepVSjet_2);
   fChain->SetBranchAddress("byVTightDeepVSjet_2", &byVTightDeepVSjet_2, &b_byVTightDeepVSjet_2);
   fChain->SetBranchAddress("byVVTightDeepVSjet_2", &byVVTightDeepVSjet_2, &b_byVVTightDeepVSjet_2);
   fChain->SetBranchAddress("byVVVLooseDeepVSmu_2", &byVVVLooseDeepVSmu_2, &b_byVVVLooseDeepVSmu_2);
   fChain->SetBranchAddress("byVVLooseDeepVSmu_2", &byVVLooseDeepVSmu_2, &b_byVVLooseDeepVSmu_2);
   fChain->SetBranchAddress("byVLooseDeepVSmu_2", &byVLooseDeepVSmu_2, &b_byVLooseDeepVSmu_2);
   fChain->SetBranchAddress("byLooseDeepVSmu_2", &byLooseDeepVSmu_2, &b_byLooseDeepVSmu_2);
   fChain->SetBranchAddress("byMediumDeepVSmu_2", &byMediumDeepVSmu_2, &b_byMediumDeepVSmu_2);
   fChain->SetBranchAddress("byTightDeepVSmu_2", &byTightDeepVSmu_2, &b_byTightDeepVSmu_2);
   fChain->SetBranchAddress("byVTightDeepVSmu_2", &byVTightDeepVSmu_2, &b_byVTightDeepVSmu_2);
   fChain->SetBranchAddress("byVVTightDeepVSmu_2", &byVVTightDeepVSmu_2, &b_byVVTightDeepVSmu_2);
   fChain->SetBranchAddress("byVVVLooseDeepVSe_2", &byVVVLooseDeepVSe_2, &b_byVVVLooseDeepVSe_2);
   fChain->SetBranchAddress("byVVLooseDeepVSe_2", &byVVLooseDeepVSe_2, &b_byVVLooseDeepVSe_2);
   fChain->SetBranchAddress("byVLooseDeepVSe_2", &byVLooseDeepVSe_2, &b_byVLooseDeepVSe_2);
   fChain->SetBranchAddress("byLooseDeepVSe_2", &byLooseDeepVSe_2, &b_byLooseDeepVSe_2);
   fChain->SetBranchAddress("byMediumDeepVSe_2", &byMediumDeepVSe_2, &b_byMediumDeepVSe_2);
   fChain->SetBranchAddress("byTightDeepVSe_2", &byTightDeepVSe_2, &b_byTightDeepVSe_2);
   fChain->SetBranchAddress("byVTightDeepVSe_2", &byVTightDeepVSe_2, &b_byVTightDeepVSe_2);
   fChain->SetBranchAddress("byVVTightDeepVSe_2", &byVVTightDeepVSe_2, &b_byVVTightDeepVSe_2);
   fChain->SetBranchAddress("numGenJets", &numGenJets, &b_numGenJets);
   fChain->SetBranchAddress("jetPt_2", &jetPt_2, &b_jetPt_2);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibReducedMINIAODFilter", &Flag_ecalBadCalibReducedMINIAODFilter, &b_Flag_ecalBadCalibReducedMINIAODFilter);
   fChain->SetBranchAddress("Flag_badMuons", &Flag_badMuons, &b_Flag_badMuons);
   fChain->SetBranchAddress("Flag_duplicateMuons", &Flag_duplicateMuons, &b_Flag_duplicateMuons);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("matchEmbFilter_Ele24Tau30_2", &matchEmbFilter_Ele24Tau30_2, &b_matchEmbFilter_Ele24Tau30_2);
   fChain->SetBranchAddress("matchEmbFilter_Ele27_1", &matchEmbFilter_Ele27_1, &b_matchEmbFilter_Ele27_1);
   fChain->SetBranchAddress("matchEmbFilter_Ele32DoubleL1v2_1", &matchEmbFilter_Ele32DoubleL1v2_1, &b_matchEmbFilter_Ele32DoubleL1v2_1);
   fChain->SetBranchAddress("matchEmbFilter_Ele32DoubleL1v1_1", &matchEmbFilter_Ele32DoubleL1v1_1, &b_matchEmbFilter_Ele32DoubleL1v1_1);
   fChain->SetBranchAddress("matchEmbFilter_Ele32_1", &matchEmbFilter_Ele32_1, &b_matchEmbFilter_Ele32_1);
   fChain->SetBranchAddress("matchEmbFilter_Ele35_1", &matchEmbFilter_Ele35_1, &b_matchEmbFilter_Ele35_1);
   fChain->SetBranchAddress("matchEmbFilter_Ele24Tau30_1", &matchEmbFilter_Ele24Tau30_1, &b_matchEmbFilter_Ele24Tau30_1);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metSig", &metSig, &b_metSig);
   fChain->SetBranchAddress("metcov00", &metcov00, &b_metcov00);
   fChain->SetBranchAddress("metcov10", &metcov10, &b_metcov10);
   fChain->SetBranchAddress("metcov11", &metcov11, &b_metcov11);
   fChain->SetBranchAddress("metcov01", &metcov01, &b_metcov01);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("met_py", &met_py, &b_met_py);
   fChain->SetBranchAddress("met_px", &met_px, &b_met_px);
   fChain->SetBranchAddress("met_UESUp", &met_UESUp, &b_met_UESUp);
   fChain->SetBranchAddress("metphi_UESUp", &metphi_UESUp, &b_metphi_UESUp);
   fChain->SetBranchAddress("met_UESDown", &met_UESDown, &b_met_UESDown);
   fChain->SetBranchAddress("metphi_UESDown", &metphi_UESDown, &b_metphi_UESDown);
   fChain->SetBranchAddress("met_responseUp", &met_responseUp, &b_met_responseUp);
   fChain->SetBranchAddress("met_responseDown", &met_responseDown, &b_met_responseDown);
   fChain->SetBranchAddress("met_resolutionUp", &met_resolutionUp, &b_met_resolutionUp);
   fChain->SetBranchAddress("met_resolutionDown", &met_resolutionDown, &b_met_resolutionDown);
   fChain->SetBranchAddress("metphi_responseUp", &metphi_responseUp, &b_metphi_responseUp);
   fChain->SetBranchAddress("metphi_responseDown", &metphi_responseDown, &b_metphi_responseDown);
   fChain->SetBranchAddress("metphi_resolutionUp", &metphi_resolutionUp, &b_metphi_resolutionUp);
   fChain->SetBranchAddress("metphi_resolutionDown", &metphi_resolutionDown, &b_metphi_resolutionDown);
   fChain->SetBranchAddress("met_JetEta0to3Up", &met_JetEta0to3Up, &b_met_JetEta0to3Up);
   fChain->SetBranchAddress("metphi_JetEta0to3Up", &metphi_JetEta0to3Up, &b_metphi_JetEta0to3Up);
   fChain->SetBranchAddress("met_JetEta0to3Down", &met_JetEta0to3Down, &b_met_JetEta0to3Down);
   fChain->SetBranchAddress("metphi_JetEta0to3Down", &metphi_JetEta0to3Down, &b_metphi_JetEta0to3Down);
   fChain->SetBranchAddress("met_JetEC2Up", &met_JetEC2Up, &b_met_JetEC2Up);
   fChain->SetBranchAddress("metphi_JetEC2Up", &metphi_JetEC2Up, &b_metphi_JetEC2Up);
   fChain->SetBranchAddress("met_JetEC2Down", &met_JetEC2Down, &b_met_JetEC2Down);
   fChain->SetBranchAddress("metphi_JetEC2Down", &metphi_JetEC2Down, &b_metphi_JetEC2Down);
   fChain->SetBranchAddress("met_JetRelativeBalUp", &met_JetRelativeBalUp, &b_met_JetRelativeBalUp);
   fChain->SetBranchAddress("metphi_JetRelativeBalUp", &metphi_JetRelativeBalUp, &b_metphi_JetRelativeBalUp);
   fChain->SetBranchAddress("met_JetRelativeBalDown", &met_JetRelativeBalDown, &b_met_JetRelativeBalDown);
   fChain->SetBranchAddress("metphi_JetRelativeBalDown", &metphi_JetRelativeBalDown, &b_metphi_JetRelativeBalDown);
   fChain->SetBranchAddress("met_JetRelativeSampleUp", &met_JetRelativeSampleUp, &b_met_JetRelativeSampleUp);
   fChain->SetBranchAddress("metphi_JetRelativeSampleUp", &metphi_JetRelativeSampleUp, &b_metphi_JetRelativeSampleUp);
   fChain->SetBranchAddress("met_JetRelativeSampleDown", &met_JetRelativeSampleDown, &b_met_JetRelativeSampleDown);
   fChain->SetBranchAddress("metphi_JetRelativeSampleDown", &metphi_JetRelativeSampleDown, &b_metphi_JetRelativeSampleDown);
   fChain->SetBranchAddress("met_JetEta3to5Up", &met_JetEta3to5Up, &b_met_JetEta3to5Up);
   fChain->SetBranchAddress("metphi_JetEta3to5Up", &metphi_JetEta3to5Up, &b_metphi_JetEta3to5Up);
   fChain->SetBranchAddress("met_JetEta3to5Down", &met_JetEta3to5Down, &b_met_JetEta3to5Down);
   fChain->SetBranchAddress("metphi_JetEta3to5Down", &metphi_JetEta3to5Down, &b_metphi_JetEta3to5Down);
   fChain->SetBranchAddress("met_JetEta0to5Up", &met_JetEta0to5Up, &b_met_JetEta0to5Up);
   fChain->SetBranchAddress("metphi_JetEta0to5Up", &metphi_JetEta0to5Up, &b_metphi_JetEta0to5Up);
   fChain->SetBranchAddress("met_JetEta0to5Down", &met_JetEta0to5Down, &b_met_JetEta0to5Down);
   fChain->SetBranchAddress("metphi_JetEta0to5Down", &metphi_JetEta0to5Down, &b_metphi_JetEta0to5Down);
   fChain->SetBranchAddress("passEle27", &passEle27, &b_passEle27);
   fChain->SetBranchAddress("passEle32", &passEle32, &b_passEle32);
   fChain->SetBranchAddress("passEle35", &passEle35, &b_passEle35);
   fChain->SetBranchAddress("passEle24Tau30", &passEle24Tau30, &b_passEle24Tau30);
   fChain->SetBranchAddress("matchEle27_1", &matchEle27_1, &b_matchEle27_1);
   fChain->SetBranchAddress("matchEle32_1", &matchEle32_1, &b_matchEle32_1);
   fChain->SetBranchAddress("matchEle35_1", &matchEle35_1, &b_matchEle35_1);
   fChain->SetBranchAddress("matchEle24Tau30_1", &matchEle24Tau30_1, &b_matchEle24Tau30_1);
   fChain->SetBranchAddress("matchEle24Tau30_2", &matchEle24Tau30_2, &b_matchEle24Tau30_2);
   fChain->SetBranchAddress("filterEle27_1", &filterEle27_1, &b_filterEle27_1);
   fChain->SetBranchAddress("filterEle32_1", &filterEle32_1, &b_filterEle32_1);
   fChain->SetBranchAddress("filterEle35_1", &filterEle35_1, &b_filterEle35_1);
   fChain->SetBranchAddress("filterEle24Tau30_1", &filterEle24Tau30_1, &b_filterEle24Tau30_1);
   fChain->SetBranchAddress("filterEle24Tau30_2", &filterEle24Tau30_2, &b_filterEle24Tau30_2);
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("mjj_JetRelativeBalUp", &mjj_JetRelativeBalUp, &b_mjj_JetRelativeBalUp);
   fChain->SetBranchAddress("mjj_JetRelativeBalDown", &mjj_JetRelativeBalDown, &b_mjj_JetRelativeBalDown);
   fChain->SetBranchAddress("mjj_JetRelativeSampleUp", &mjj_JetRelativeSampleUp, &b_mjj_JetRelativeSampleUp);
   fChain->SetBranchAddress("mjj_JetRelativeSampleDown", &mjj_JetRelativeSampleDown, &b_mjj_JetRelativeSampleDown);
   fChain->SetBranchAddress("mjj_JetEta0to3Up", &mjj_JetEta0to3Up, &b_mjj_JetEta0to3Up);
   fChain->SetBranchAddress("mjj_JetEta0to3Down", &mjj_JetEta0to3Down, &b_mjj_JetEta0to3Down);
   fChain->SetBranchAddress("mjj_JetEC2Up", &mjj_JetEC2Up, &b_mjj_JetEC2Up);
   fChain->SetBranchAddress("mjj_JetEC2Down", &mjj_JetEC2Down, &b_mjj_JetEC2Down);
   fChain->SetBranchAddress("mjj_JetEta0to5Up", &mjj_JetEta0to5Up, &b_mjj_JetEta0to5Up);
   fChain->SetBranchAddress("mjj_JetEta0to5Down", &mjj_JetEta0to5Down, &b_mjj_JetEta0to5Down);
   fChain->SetBranchAddress("mjj_JetEta3to5Up", &mjj_JetEta3to5Up, &b_mjj_JetEta3to5Up);
   fChain->SetBranchAddress("mjj_JetEta3to5Down", &mjj_JetEta3to5Down, &b_mjj_JetEta3to5Down);
   fChain->SetBranchAddress("gen_match_1", &gen_match_1, &b_gen_match_1);
   fChain->SetBranchAddress("gen_match_2", &gen_match_2, &b_gen_match_2);
   fChain->SetBranchAddress("nbtag", &nbtag, &b_nbtag);
   fChain->SetBranchAddress("nbtagL", &nbtagL, &b_nbtagL);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njets_JetRelativeBalUp", &njets_JetRelativeBalUp, &b_njets_JetRelativeBalUp);
   fChain->SetBranchAddress("njets_JetRelativeBalDown", &njets_JetRelativeBalDown, &b_njets_JetRelativeBalDown);
   fChain->SetBranchAddress("njets_JetRelativeSampleUp", &njets_JetRelativeSampleUp, &b_njets_JetRelativeSampleUp);
   fChain->SetBranchAddress("njets_JetRelativeSampleDown", &njets_JetRelativeSampleDown, &b_njets_JetRelativeSampleDown);
   fChain->SetBranchAddress("njets_JetEta0to3Up", &njets_JetEta0to3Up, &b_njets_JetEta0to3Up);
   fChain->SetBranchAddress("njets_JetEta0to3Down", &njets_JetEta0to3Down, &b_njets_JetEta0to3Down);
   fChain->SetBranchAddress("njets_JetEC2Up", &njets_JetEC2Up, &b_njets_JetEC2Up);
   fChain->SetBranchAddress("njets_JetEC2Down", &njets_JetEC2Down, &b_njets_JetEC2Down);
   fChain->SetBranchAddress("njets_JetEta0to5Up", &njets_JetEta0to5Up, &b_njets_JetEta0to5Up);
   fChain->SetBranchAddress("njets_JetEta0to5Down", &njets_JetEta0to5Down, &b_njets_JetEta0to5Down);
   fChain->SetBranchAddress("njets_JetEta3to5Up", &njets_JetEta3to5Up, &b_njets_JetEta3to5Up);
   fChain->SetBranchAddress("njets_JetEta3to5Down", &njets_JetEta3to5Down, &b_njets_JetEta3to5Down);
   fChain->SetBranchAddress("jpt_1", &jpt_1, &b_jpt_1);
   fChain->SetBranchAddress("jeta_1", &jeta_1, &b_jeta_1);
   fChain->SetBranchAddress("jphi_1", &jphi_1, &b_jphi_1);
   fChain->SetBranchAddress("jpt_2", &jpt_2, &b_jpt_2);
   fChain->SetBranchAddress("jeta_2", &jeta_2, &b_jeta_2);
   fChain->SetBranchAddress("jphi_2", &jphi_2, &b_jphi_2);
   fChain->SetBranchAddress("jcsv_1", &jcsv_1, &b_jcsv_1);
   fChain->SetBranchAddress("jcsv_2", &jcsv_2, &b_jcsv_2);
   fChain->SetBranchAddress("jpt_JetEta0to3Up_1", &jpt_JetEta0to3Up_1, &b_jpt_JetEta0to3Up_1);
   fChain->SetBranchAddress("jpt_JetEta0to3Down_1", &jpt_JetEta0to3Down_1, &b_jpt_JetEta0to3Down_1);
   fChain->SetBranchAddress("jpt_JetEC2Up_1", &jpt_JetEC2Up_1, &b_jpt_JetEC2Up_1);
   fChain->SetBranchAddress("jpt_JetEC2Down_1", &jpt_JetEC2Down_1, &b_jpt_JetEC2Down_1);
   fChain->SetBranchAddress("jpt_JetEta3to5Up_1", &jpt_JetEta3to5Up_1, &b_jpt_JetEta3to5Up_1);
   fChain->SetBranchAddress("jpt_JetEta3to5Down_1", &jpt_JetEta3to5Down_1, &b_jpt_JetEta3to5Down_1);
   fChain->SetBranchAddress("jpt_JetEta0to5Up_1", &jpt_JetEta0to5Up_1, &b_jpt_JetEta0to5Up_1);
   fChain->SetBranchAddress("jpt_JetEta0to5Down_1", &jpt_JetEta0to5Down_1, &b_jpt_JetEta0to5Down_1);
   fChain->SetBranchAddress("jpt_JetRelativeBalUp_1", &jpt_JetRelativeBalUp_1, &b_jpt_JetRelativeBalUp_1);
   fChain->SetBranchAddress("jpt_JetRelativeBalDown_1", &jpt_JetRelativeBalDown_1, &b_jpt_JetRelativeBalDown_1);
   fChain->SetBranchAddress("jpt_JetRelativeSampleUp_1", &jpt_JetRelativeSampleUp_1, &b_jpt_JetRelativeSampleUp_1);
   fChain->SetBranchAddress("jpt_JetRelativeSampleDown_1", &jpt_JetRelativeSampleDown_1, &b_jpt_JetRelativeSampleDown_1);
   fChain->SetBranchAddress("jpt_JetEta0to3Up_2", &jpt_JetEta0to3Up_2, &b_jpt_JetEta0to3Up_2);
   fChain->SetBranchAddress("jpt_JetEta0to3Down_2", &jpt_JetEta0to3Down_2, &b_jpt_JetEta0to3Down_2);
   fChain->SetBranchAddress("jpt_JetEC2Up_2", &jpt_JetEC2Up_2, &b_jpt_JetEC2Up_2);
   fChain->SetBranchAddress("jpt_JetEC2Down_2", &jpt_JetEC2Down_2, &b_jpt_JetEC2Down_2);
   fChain->SetBranchAddress("jpt_JetEta3to5Up_2", &jpt_JetEta3to5Up_2, &b_jpt_JetEta3to5Up_2);
   fChain->SetBranchAddress("jpt_JetEta3to5Down_2", &jpt_JetEta3to5Down_2, &b_jpt_JetEta3to5Down_2);
   fChain->SetBranchAddress("jpt_JetEta0to5Up_2", &jpt_JetEta0to5Up_2, &b_jpt_JetEta0to5Up_2);
   fChain->SetBranchAddress("jpt_JetEta0to5Down_2", &jpt_JetEta0to5Down_2, &b_jpt_JetEta0to5Down_2);
   fChain->SetBranchAddress("jpt_JetRelativeBalUp_2", &jpt_JetRelativeBalUp_2, &b_jpt_JetRelativeBalUp_2);
   fChain->SetBranchAddress("jpt_JetRelativeBalDown_2", &jpt_JetRelativeBalDown_2, &b_jpt_JetRelativeBalDown_2);
   fChain->SetBranchAddress("jpt_JetRelativeSampleUp_2", &jpt_JetRelativeSampleUp_2, &b_jpt_JetRelativeSampleUp_2);
   fChain->SetBranchAddress("jpt_JetRelativeSampleDown_2", &jpt_JetRelativeSampleDown_2, &b_jpt_JetRelativeSampleDown_2);
   fChain->SetBranchAddress("bpt_1", &bpt_1, &b_bpt_1);
   fChain->SetBranchAddress("bflavor_1", &bflavor_1, &b_bflavor_1);
   fChain->SetBranchAddress("beta_1", &beta_1, &b_beta_1);
   fChain->SetBranchAddress("bphi_1", &bphi_1, &b_bphi_1);
   fChain->SetBranchAddress("bpt_2", &bpt_2, &b_bpt_2);
   fChain->SetBranchAddress("bflavor_2", &bflavor_2, &b_bflavor_2);
   fChain->SetBranchAddress("beta_2", &beta_2, &b_beta_2);
   fChain->SetBranchAddress("bphi_2", &bphi_2, &b_bphi_2);
   fChain->SetBranchAddress("pt_top1", &pt_top1, &b_pt_top1);
   fChain->SetBranchAddress("pt_top2", &pt_top2, &b_pt_top2);
   fChain->SetBranchAddress("genweight", &genweight, &b_genweight);
   fChain->SetBranchAddress("gen_Higgs_pt", &gen_Higgs_pt, &b_gen_Higgs_pt);
   fChain->SetBranchAddress("gen_Higgs_mass", &gen_Higgs_mass, &b_gen_Higgs_mass);
   Notify();
}

Bool_t smhet_2017::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void smhet_2017::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t smhet_2017::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef smhet_2017_cxx
