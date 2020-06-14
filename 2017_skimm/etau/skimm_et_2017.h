//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 14 01:21:53 2018 by ROOT version 6.10/09
// from TTree EventTree/Event data (tag V08_00_26_03)
// found on file: /hdfs/store/user/jmadhusu/MonoHiggs_MC2017/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50/180603_140815/0000/ggtree_mc_1.root
//////////////////////////////////////////////////////////

#ifndef skimm_et_2017_h
#define skimm_et_2017_h

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
//#include <TDCacheFile.h>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "TString.h"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"


using namespace std;

class skimm_et_2017 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   //TTree *outTree=new TTree("eventTree","eventTree");
   TFile *fileName;
   TTree *tree;
   TTree *newtree;
   TH1F* h_nEvents;
   

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   Float_t         vtxX;
   Float_t         vtxY;
   Float_t         vtxZ;
   Int_t           vtxNtrks;
   Bool_t          vtx_isFake;
   Int_t           vtx_ndof;
   Float_t         vtx_rho;
   Bool_t          isGoodVtx;
   Int_t           nGoodVtx;
   Float_t         rho;
   Float_t         rhoCentral;
   Double_t        prefiringweight;
   Double_t        prefiringweightup;
   Double_t        prefiringweightdown;
   ULong64_t       HLTEleMuX;
   ULong64_t       HLTEleMuXIsPrescaled;
   ULong64_t       HLTEleMuXRejectedByPS;
   ULong64_t       HLTPho;
   ULong64_t       HLTPhoIsPrescaled;
   ULong64_t       HLTPhoRejectedByPS;
   ULong64_t       HLTTau;
   ULong64_t       HLTTauIsPrescaled;
   ULong64_t       HLTTauRejectedByPS;
   ULong64_t       HLTMet;
   ULong64_t       HLTMetIsPrescaled;
   ULong64_t       HLTMetRejectedByPS;
   ULong64_t       HLTJet;
   ULong64_t       HLTJetIsPrescaled;
   ULong64_t       HLTJetRejectedByPS;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoUnCalibE;
   vector<float>   *phoUnCalibESigma;
   vector<float>   *phoCalibE;
   vector<float>   *phoCalibESigma;
   vector<float>   *phoCalibEt;
   vector<float>   *phoEnergyScale;
   vector<float>   *phoEnergySigma;
   vector<float>   *phoSCE;
   vector<float>   *phoSCRawE;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<int>     *phohasPixelSeed;
   vector<int>     *phoEleVeto;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoHoverE;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoSigmaIEtaIPhiFull5x5;
   vector<float>   *phoSigmaIPhiIPhiFull5x5;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFChWorstIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoIDMVA;
   vector<unsigned short> *phoIDbit;
   vector<float>   *phoSeedTime;
   vector<float>   *phoSeedEnergy;
   vector<ULong64_t> *phoFiredSingleTrgs;
   vector<ULong64_t> *phoFiredDoubleTrgs;
   vector<ULong64_t> *phoFiredTripleTrgs;
   vector<ULong64_t> *phoFiredL1Trgs;
   vector<float>   *phoScale_up;
   vector<float>   *phoScale_dn;
   vector<float>   *phoScale_stat_up;
   vector<float>   *phoScale_stat_dn;
   vector<float>   *phoScale_syst_up;
   vector<float>   *phoScale_syst_dn;
   vector<float>   *phoScale_gain_up;
   vector<float>   *phoScale_gain_dn;
   vector<float>   *phoResol_up;
   vector<float>   *phoResol_dn;
   vector<float>   *phoResol_rho_up;
   vector<float>   *phoResol_rho_dn;
   vector<float>   *phoResol_phi_up;
   vector<float>   *phoResol_phi_dn;
   Int_t           nJet;
   vector<float>   *jetPt;
   vector<float>   *jetE;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetRawPt;
   vector<float>   *jetRawE;
   vector<float>   *jetMt;
   vector<float>   *jetArea;
   vector<float>   *jetMass;
   vector<float>   *jetMaxDistance;
   vector<float>   *jetPhiPhiMoment;
   vector<float>   *jetConstituentEtaPhiSpread;
   vector<float>   *jetConstituentPtDistribution;
   vector<float>   *jetPileup;
   vector<unsigned short> *jetID;
   vector<float>   *jetPUID;
   vector<int>     *jetPUFullID;
   vector<int>     *jetPartonID;
   vector<int>     *jetHadFlvr;
   vector<float>   *jetJECUnc;
   vector<float>   *jetCEF;
   vector<float>   *jetNEF;
   vector<float>   *jetCHF;
   vector<float>   *jetNHF;
   vector<int>     *jetPhotonEnF;
   vector<int>     *jetElectronEnF;
   vector<float>   *jetMuonEnF;
   vector<float>   *jetChargedMuonEnF;
   vector<float>   *jetHFHAE;
   vector<float>   *jetHFEME;
   vector<int>     *jetNConst;
   vector<int>     *jetNConstituents;
   vector<int>     *jetNCharged;
   vector<int>     *jetNNeutral;
   vector<int>     *jetNChargedHad;
   vector<int>     *jetNNeutralHad;
   vector<int>     *jetNPhoton;
   vector<int>     *jetNElectron;
   vector<int>     *jetNMuon;
   vector<float>   *jetCSV2BJetTags;
   vector<float>   *jetDeepCSVTags_b;
   vector<float>   *jetDeepCSVTags_bb;
   vector<float>   *jetDeepCSVTags_c;
   vector<float>   *jetDeepCSVTags_udsg;
   vector<float>   *jetDeepFlavour_b;
   vector<float>   *jetDeepFlavour_bb;
   vector<float>   *jetDeepFlavour_lepb;
   vector<float>   *jetDeepFlavour_c;
   vector<float>   *jetDeepFlavour_uds;
   vector<float>   *jetDeepFlavour_g;
   vector<float>   *jetetaWidth;
   vector<float>   *jetphiWidth;
   vector<vector<float> > *jetConstPt;
   vector<vector<float> > *jetConstEt;
   vector<vector<float> > *jetConstEta;
   vector<vector<float> > *jetConstPhi;
   vector<vector<int> > *jetConstPdgId;
   vector<float>   *jetGenJetE;
   vector<float>   *jetGenJetPt;
   vector<float>   *jetGenJetEta;
   vector<float>   *jetGenJetPhi;
   vector<int>     *jetGenPartonID;
   vector<float>   *jetGenE;
   vector<float>   *jetGenPt;
   vector<float>   *jetGenEta;
   vector<float>   *jetGenPhi;
   vector<int>     *jetGenPartonMomID;
   vector<float>   *jetP4Smear;
   vector<float>   *jetP4SmearUp;
   vector<float>   *jetP4SmearDo;
   Int_t           nEle;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9Full5x5;
   vector<float>   *eleE;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleSIP;
   vector<float>   *eleUnCalibE;
   vector<float>   *eleUnCalibESigma;
   vector<float>   *eleCalibEecalonly;
   vector<float>   *eleCalibE;
   vector<float>   *eleCalibESigma;
   vector<float>   *eleCalibEt;
   vector<float>   *eleCalibEtSigma;
   vector<float>   *eleEnergyScale;
   vector<float>   *eleEnergySigma;
   vector<float>   *eleSCRawE;
   vector<float>   *eleSCE;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCEtaWidth;
   vector<float>   *eleSCPhiWidth;
   vector<float>   *eleHoverE;
   vector<float>   *eleEoverP;
   vector<float>   *eleEoverPout;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eledEtaAtCalo;
   vector<float>   *eledEtaseedAtVtx;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleSigmaIPhiIPhiFull5x5;
   vector<int>     *eleConvVeto;
   vector<int>     *eleMissHits;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *elePFClusEcalIso;
   vector<float>   *elePFClusHcalIso;
   vector<ULong64_t> *eleFiredSingleTrgs;
   vector<ULong64_t> *eleFiredDoubleTrgs;
   vector<ULong64_t> *eleFiredL1Trgs;
   vector<float>   *eleHEEPID;
   vector<float>   *eleMVAIsoID;
   vector<float>   *eleMVAnoIsoID;
   vector<unsigned short> *eleIDbit;
   vector<float>   *eleTrkdxy;
   vector<float>   *eleKFHits;
   vector<float>   *eleKFChi2;
   vector<float>   *eleGSFChi2;
   vector<float>   *eleScale_up;
   vector<float>   *eleScale_dn;
   vector<float>   *eleScale_stat_up;
   vector<float>   *eleScale_stat_dn;
   vector<float>   *eleScale_syst_up;
   vector<float>   *eleScale_syst_dn;
   vector<float>   *eleScale_gain_up;
   vector<float>   *eleScale_gain_dn;
   vector<float>   *eleResol_up;
   vector<float>   *eleResol_dn;
   vector<float>   *eleResol_rho_up;
   vector<float>   *eleResol_rho_dn;
   vector<float>   *eleResol_phi_up;
   vector<float>   *eleResol_phi_dn;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muE;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<int>     *muCharge;
   vector<int>     *muType;
   vector<unsigned short> *muIDbit;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muSIP;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelLayers;
   vector<int>     *muPixelHits;
   vector<int>     *muMuonHits;
   vector<int>     *muStations;
   vector<int>     *muMatches;
   vector<int>     *muTrkQuality;
   vector<float>   *muInnervalidFraction;
   vector<float>   *muIsoTrk;
   vector<float>   *muPFChIso;
   vector<float>   *muPFPhoIso;
   vector<float>   *muPFNeuIso;
   vector<float>   *muPFPUIso;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muBestTrkPtError;
   vector<float>   *muBestTrkPt;
   vector<int>     *muBestTrkType;
   vector<ULong64_t> *muFiredTrgs;
   vector<ULong64_t> *muFiredL1Trgs;
   Int_t           nTau;
   vector<float>   *tau_Pt;
   vector<float>   *tau_Et;
   vector<float>   *tau_Eta;
   vector<float>   *tau_Phi;
   vector<float>   *tau_Charge;
   vector<int>     *tau_DecayMode;
   vector<int>     *tau_decayModeFindingNewDMs;
   vector<float>   *tau_P;
   vector<float>   *tau_Vz;
   vector<float>   *tau_Energy;
   vector<float>   *tau_Mass;
   vector<float>   *tau_Dxy;
   vector<float>   *tau_ZImpact;
   vector<float>   *tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<float>   *tau_chargedIsoPtSum;
   vector<float>   *tau_neutralIsoPtSum;
   vector<float>   *tau_neutralIsoPtSumWeight;
   vector<float>   *tau_footprintCorrection;
   vector<float>   *tau_photonPtSumOutsideSignalCone;
   vector<float>   *tau_puCorrPtSum;
   vector<int>     *tau_NumSignalPFChargedHadrCands;
   vector<int>     *tau_NumSignalPFNeutrHadrCands;
   vector<int>     *tau_NumSignalPFGammaCands;
   vector<int>     *tau_NumSignalPFCands;
   vector<int>     *tau_NumIsolationPFChargedHadrCands;
   vector<int>     *tau_NumIsolationPFNeutrHadrCands;
   vector<int>     *tau_NumIsolationPFGammaCands;
   vector<int>     *tau_NumIsolationPFCands;
   vector<float>   *tau_LeadChargedHadronEta;
   vector<float>   *tau_LeadChargedHadronPhi;
   vector<float>   *tau_LeadChargedHadronPt;
   vector<float>   *tau_LeadChargedHadron_dz;
   vector<float>   *tau_LeadChargedHadron_dxy;
   vector<unsigned int> *tau_IDbits;
   vector<float>   *tau_byIsolationMVArun2017v2DBoldDMwLTraw2017;
   vector<float>   *tau_byVVVLooseDeepTau2017v2p1VSjet;
   vector<float>   *tau_byVVLooseDeepTau2017v2p1VSjet;
   vector<float>   *tau_byVLooseDeepTau2017v2p1VSjet;
   vector<float>   *tau_byLooseDeepTau2017v2p1VSjet;
   vector<float>   *tau_byMediumDeepTau2017v2p1VSjet;
   vector<float>   *tau_byTightDeepTau2017v2p1VSjet;
   vector<float>   *tau_byVTightDeepTau2017v2p1VSjet;
   vector<float>   *tau_byVVTightDeepTau2017v2p1VSjet;
   vector<float>   *tau_byVVVLooseDeepTau2017v2p1VSe;
   vector<float>   *tau_byVVLooseDeepTau2017v2p1VSe;
   vector<float>   *tau_byVLooseDeepTau2017v2p1VSe;
   vector<float>   *tau_byLooseDeepTau2017v2p1VSe;
   vector<float>   *tau_byMediumDeepTau2017v2p1VSe;
   vector<float>   *tau_byTightDeepTau2017v2p1VSe;
   vector<float>   *tau_byVTightDeepTau2017v2p1VSe;
   vector<float>   *tau_byVVTightDeepTau2017v2p1VSe;
   vector<float>   *tau_byVLooseDeepTau2017v2p1VSmu;
   vector<float>   *tau_byLooseDeepTau2017v2p1VSmu;
   vector<float>   *tau_byMediumDeepTau2017v2p1VSmu;
   vector<float>   *tau_byTightDeepTau2017v2p1VSmu;
   vector<ULong64_t> *tauFiredTrgs;
   vector<ULong64_t> *tauFiredL1Trgs;
   Float_t         genMET;
   Float_t         genMETPhi;
   UShort_t        metFilters;
   Float_t         caloMET;
   Float_t         caloMETPhi;
   Float_t         caloMETsumEt;
   Float_t         pfMETCorr;
   Float_t         pfMETPhiCorr;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         pfMET_T1JERUp;
   Float_t         pfMET_T1JERDo;
   Float_t         pfMET_T1JESUp;
   Float_t         pfMET_T1JESDo;
   Float_t         pfMET_T1UESUp;
   Float_t         pfMET_T1UESDo;
   Float_t         pfMETPhi_T1JESUp;
   Float_t         pfMETPhi_T1JESDo;
   Float_t         pfMETPhi_T1UESUp;
   Float_t         pfMETPhi_T1UESDo;
   vector<float>   *pdf;
   Float_t         pthat;
   Float_t         processID;
   Float_t         genWeight;
   Float_t         genHT;
   Float_t         pdfWeight;
   vector<float>   *pdfSystWeight;
   Int_t           nPUInfo;
   vector<int>     *nPU;
   vector<int>     *puBX;
   vector<float>   *puTrue;
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<float>   *mcVtx;
   vector<float>   *mcVty;
   vector<float>   *mcVtz;
   vector<float>   *mcPt;
   vector<float>   *mcMass;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<int>     *mcStatus;
   vector<unsigned short> *mcStatusFlag;
   vector<int>     *mcIndex;
   vector<int>     *mcDaughterPID;
   vector<float>   *mcCharge;
   vector<int>     *mcMotherPID;
   vector<int>     *mcMotherIndex;
   vector<int>     *mcMotherStatus;
   vector<int>     *mcDaughterStatus;
   vector<int>     *mcDaughterList;
   vector<unsigned short> *mcTauDecayMode;
   vector<unsigned short> *genMatch2;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtxX;   //!
   TBranch        *b_vtxY;   //!
   TBranch        *b_vtxZ;   //!
   TBranch        *b_vtxNtrks;   //!
   TBranch        *b_vtx_isFake;   //!
   TBranch        *b_vtx_ndof;   //!
   TBranch        *b_vtx_rho;   //!
   TBranch        *b_isGoodVtx;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_prefiringweight;   //!
   TBranch        *b_prefiringweightup;   //!
   TBranch        *b_prefiringweightdown;   //!
   TBranch        *b_HLTEleMuX;   //!
   TBranch        *b_HLTEleMuXIsPrescaled;   //!
   TBranch        *b_HLTEleMuXRejectedByPS;   //!
   TBranch        *b_HLTPho;   //!
   TBranch        *b_HLTPhoIsPrescaled;   //!
   TBranch        *b_HLTPhoRejectedByPS;   //!
   TBranch        *b_HLTTau;   //!
   TBranch        *b_HLTTauIsPrescaled;   //!
   TBranch        *b_HLTTauRejectedByPS;   //!
   TBranch        *b_HLTMet;   //!
   TBranch        *b_HLTMetIsPrescaled;   //!
   TBranch        *b_HLTMetRejectedByPS;   //!
   TBranch        *b_HLTJet;   //!
   TBranch        *b_HLTJetIsPrescaled;   //!
   TBranch        *b_HLTJetRejectedByPS;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoUnCalibE;   //!
   TBranch        *b_phoUnCalibESigma;   //!
   TBranch        *b_phoCalibE;   //!
   TBranch        *b_phoCalibESigma;   //!
   TBranch        *b_phoCalibEt;   //!
   TBranch        *b_phoEnergyScale;   //!
   TBranch        *b_phoEnergySigma;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
   TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFChWorstIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoIDbit;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoFiredSingleTrgs;   //!
   TBranch        *b_phoFiredDoubleTrgs;   //!
   TBranch        *b_phoFiredTripleTrgs;   //!
   TBranch        *b_phoFiredL1Trgs;   //!
   TBranch        *b_phoScale_up;   //!
   TBranch        *b_phoScale_dn;   //!
   TBranch        *b_phoScale_stat_up;   //!
   TBranch        *b_phoScale_stat_dn;   //!
   TBranch        *b_phoScale_syst_up;   //!
   TBranch        *b_phoScale_syst_dn;   //!
   TBranch        *b_phoScale_gain_up;   //!
   TBranch        *b_phoScale_gain_dn;   //!
   TBranch        *b_phoResol_up;   //!
   TBranch        *b_phoResol_dn;   //!
   TBranch        *b_phoResol_rho_up;   //!
   TBranch        *b_phoResol_rho_dn;   //!
   TBranch        *b_phoResol_phi_up;   //!
   TBranch        *b_phoResol_phi_dn;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawE;   //!
   TBranch        *b_jetMt;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetMaxDistance;   //!
   TBranch        *b_jetPhiPhiMoment;   //!
   TBranch        *b_jetConstituentEtaPhiSpread;   //!
   TBranch        *b_jetConstituentPtDistribution;   //!
   TBranch        *b_jetPileup;   //!
   TBranch        *b_jetID;   //!
   TBranch        *b_jetPUID;   //!
   TBranch        *b_jetPUFullID;   //!
   TBranch        *b_jetPartonID;   //!
   TBranch        *b_jetHadFlvr;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetPhotonEnF;   //!
   TBranch        *b_jetElectronEnF;   //!
   TBranch        *b_jetMuonEnF;   //!
   TBranch        *b_jetChargedMuonEnF;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConst;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetNCharged;   //!
   TBranch        *b_jetNNeutral;   //!
   TBranch        *b_jetNChargedHad;   //!
   TBranch        *b_jetNNeutralHad;   //!
   TBranch        *b_jetNPhoton;   //!
   TBranch        *b_jetNElectron;   //!
   TBranch        *b_jetNMuon;   //!
   TBranch        *b_jetCSV2BJetTags;   //!
   TBranch        *b_jetDeepCSVTags_b;   //!
   TBranch        *b_jetDeepCSVTags_bb;   //!
   TBranch        *b_jetDeepCSVTags_c;   //!
   TBranch        *b_jetDeepCSVTags_udsg;   //!
   TBranch        *b_jetDeepFlavour_b;   //!
   TBranch        *b_jetDeepFlavour_bb;   //!
   TBranch        *b_jetDeepFlavour_lepb;   //!
   TBranch        *b_jetDeepFlavour_c;   //!
   TBranch        *b_jetDeepFlavour_uds;   //!
   TBranch        *b_jetDeepFlavour_g;   //!
   TBranch        *b_jetetaWidth;   //!
   TBranch        *b_jetphiWidth;   //!
   TBranch        *b_jetConstPt;   //!
   TBranch        *b_jetConstEt;   //!
   TBranch        *b_jetConstEta;   //!
   TBranch        *b_jetConstPhi;   //!
   TBranch        *b_jetConstPdgId;   //!
   TBranch        *b_jetGenJetE;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenE;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_jetGenPartonMomID;   //!
   TBranch        *b_jetP4Smear;   //!
   TBranch        *b_jetP4SmearUp;   //!
   TBranch        *b_jetP4SmearDo;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleE;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleSIP;   //!
   TBranch        *b_eleUnCalibE;   //!
   TBranch        *b_eleUnCalibESigma;   //!
   TBranch        *b_eleCalibEecalonly;   //!
   TBranch        *b_eleCalibE;   //!
   TBranch        *b_eleCalibESigma;   //!
   TBranch        *b_eleCalibEt;   //!
   TBranch        *b_eleCalibEtSigma;   //!
   TBranch        *b_eleEnergyScale;   //!
   TBranch        *b_eleEnergySigma;   //!
   TBranch        *b_eleSCRawE;   //!
   TBranch        *b_eleSCE;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPout;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eledEtaAtCalo;   //!
   TBranch        *b_eledEtaseedAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_elePFClusEcalIso;   //!
   TBranch        *b_elePFClusHcalIso;   //!
   TBranch        *b_eleFiredSingleTrgs;   //!
   TBranch        *b_eleFiredDoubleTrgs;   //!
   TBranch        *b_eleFiredL1Trgs;   //!
   TBranch        *b_eleHEEPID;   //!
   TBranch        *b_eleMVAIsoID;   //!
   TBranch        *b_eleMVAnoIsoID;   //!
   TBranch        *b_eleIDbit;   //!
   TBranch        *b_eleTrkdxy;   //!
   TBranch        *b_eleKFHits;   //!
   TBranch        *b_eleKFChi2;   //!
   TBranch        *b_eleGSFChi2;   //!
   TBranch        *b_eleScale_up;   //!
   TBranch        *b_eleScale_dn;   //!
   TBranch        *b_eleScale_stat_up;   //!
   TBranch        *b_eleScale_stat_dn;   //!
   TBranch        *b_eleScale_syst_up;   //!
   TBranch        *b_eleScale_syst_dn;   //!
   TBranch        *b_eleScale_gain_up;   //!
   TBranch        *b_eleScale_gain_dn;   //!
   TBranch        *b_eleResol_up;   //!
   TBranch        *b_eleResol_dn;   //!
   TBranch        *b_eleResol_rho_up;   //!
   TBranch        *b_eleResol_rho_dn;   //!
   TBranch        *b_eleResol_phi_up;   //!
   TBranch        *b_eleResol_phi_dn;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muE;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIDbit;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muSIP;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muMatches;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muBestTrkPtError;   //!
   TBranch        *b_muBestTrkPt;   //!
   TBranch        *b_muBestTrkType;   //!
   TBranch        *b_muFiredTrgs;   //!
   TBranch        *b_muFiredL1Trgs;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_tau_Pt;   //!
   TBranch        *b_tau_Et;   //!
   TBranch        *b_tau_Eta;   //!
   TBranch        *b_tau_Phi;   //!
   TBranch        *b_tau_Charge;   //!
   TBranch        *b_tau_DecayMode;   //!
   TBranch        *b_tau_decayModeFindingNewDMs;   //!
   TBranch        *b_tau_P;   //!
   TBranch        *b_tau_Vz;   //!
   TBranch        *b_tau_Energy;   //!
   TBranch        *b_tau_Mass;   //!
   TBranch        *b_tau_Dxy;   //!
   TBranch        *b_tau_ZImpact;   //!
   TBranch        *b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_tau_chargedIsoPtSum;   //!
   TBranch        *b_tau_neutralIsoPtSum;   //!
   TBranch        *b_tau_neutralIsoPtSumWeight;   //!
   TBranch        *b_tau_footprintCorrection;   //!
   TBranch        *b_tau_photonPtSumOutsideSignalCone;   //!
   TBranch        *b_tau_puCorrPtSum;   //!
   TBranch        *b_tau_NumSignalPFChargedHadrCands;   //!
   TBranch        *b_tau_NumSignalPFNeutrHadrCands;   //!
   TBranch        *b_tau_NumSignalPFGammaCands;   //!
   TBranch        *b_tau_NumSignalPFCands;   //!
   TBranch        *b_tau_NumIsolationPFChargedHadrCands;   //!
   TBranch        *b_tau_NumIsolationPFNeutrHadrCands;   //!
   TBranch        *b_tau_NumIsolationPFGammaCands;   //!
   TBranch        *b_tau_NumIsolationPFCands;   //!
   TBranch        *b_tau_LeadChargedHadronEta;   //!
   TBranch        *b_tau_LeadChargedHadronPhi;   //!
   TBranch        *b_tau_LeadChargedHadronPt;   //!
   TBranch        *b_tau_LeadChargedHadron_dz;   //!
   TBranch        *b_tau_LeadChargedHadron_dxy;   //!
   TBranch        *b_tau_IDbits;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBoldDMwLTraw2017;   //!
   TBranch        *b_tau_byVVVLooseDeepTau2017v2p1VSjet;   //!
   TBranch        *b_tau_byVVLooseDeepTau2017v2p1VSjet;   //!
   TBranch        *b_tau_byVLooseDeepTau2017v2p1VSjet;   //!
   TBranch        *b_tau_byLooseDeepTau2017v2p1VSjet;   //!
   TBranch        *b_tau_byMediumDeepTau2017v2p1VSjet;   //!
   TBranch        *b_tau_byTightDeepTau2017v2p1VSjet;   //!
   TBranch        *b_tau_byVTightDeepTau2017v2p1VSjet;   //!
   TBranch        *b_tau_byVVTightDeepTau2017v2p1VSjet;   //!
   TBranch        *b_tau_byVVVLooseDeepTau2017v2p1VSe;   //!
   TBranch        *b_tau_byVVLooseDeepTau2017v2p1VSe;   //!
   TBranch        *b_tau_byVLooseDeepTau2017v2p1VSe;   //!
   TBranch        *b_tau_byLooseDeepTau2017v2p1VSe;   //!
   TBranch        *b_tau_byMediumDeepTau2017v2p1VSe;   //!
   TBranch        *b_tau_byTightDeepTau2017v2p1VSe;   //!
   TBranch        *b_tau_byVTightDeepTau2017v2p1VSe;   //!
   TBranch        *b_tau_byVVTightDeepTau2017v2p1VSe;   //!
   TBranch        *b_tau_byVLooseDeepTau2017v2p1VSmu;   //!
   TBranch        *b_tau_byLooseDeepTau2017v2p1VSmu;   //!
   TBranch        *b_tau_byMediumDeepTau2017v2p1VSmu;   //!
   TBranch        *b_tau_byTightDeepTau2017v2p1VSmu;   //!
   TBranch        *b_tauFiredTrgs;   //!
   TBranch        *b_tauFiredL1Trgs;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_caloMET;   //!
   TBranch        *b_caloMETPhi;   //!
   TBranch        *b_caloMETsumEt;   //!
   TBranch        *b_pfMETCorr;   //!
   TBranch        *b_pfMETPhiCorr;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_pfMET_T1JERUp;   //!
   TBranch        *b_pfMET_T1JERDo;   //!
   TBranch        *b_pfMET_T1JESUp;   //!
   TBranch        *b_pfMET_T1JESDo;   //!
   TBranch        *b_pfMET_T1UESUp;   //!
   TBranch        *b_pfMET_T1UESDo;   //!
   TBranch        *b_pfMETPhi_T1JESUp;   //!
   TBranch        *b_pfMETPhi_T1JESDo;   //!
   TBranch        *b_pfMETPhi_T1UESUp;   //!
   TBranch        *b_pfMETPhi_T1UESDo;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genHT;   //!
   TBranch        *b_pdfWeight;   //!
   TBranch        *b_pdfSystWeight;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcVty;   //!
   TBranch        *b_mcVtz;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcStatusFlag;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcDaughterPID;   //!
   TBranch        *b_mcCharge;   //!
   TBranch        *b_mcMotherPID;   //!
   TBranch        *b_mcMotherIndex;   //!
   TBranch        *b_mcMotherStatus;   //!
   TBranch        *b_mcDaughterStatus;   //!
   TBranch        *b_mcDaughterList;   //!
   TBranch        *b_mcTauDecayMode;   //!
   TBranch        *b_genMatch2;   //!

   //   skimm_et_2017(TTree *tree=0);
   skimm_et_2017(const char* file1, const char* file2, string isMC);
   virtual ~skimm_et_2017();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree, string isMC);
   //virtual void     Init(TTree *tree);
   virtual void     Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void BookHistos(const char* file2);
   virtual void fillHistos(int histoNumber, double event_weight,int higgs_index);
   virtual bool skimming_Htt();
   virtual vector<int> skimmed_Ele(); 
   virtual vector<int> skimmed_Tau();
   virtual void fillOutTree();
   virtual double DeltaPhi(double phi1, double phi2);
   virtual double dR(int mu_index, int tau_index);

};

#endif

#ifdef skimm_et_2017_cxx
skimm_et_2017::skimm_et_2017(const char* file1, const char* file2, string isMC)
{
  TChain *chain = new TChain("phoJetNtuplizer/eventTree");
  //TH1F *inspected_events = new TH1F("inspected_events","inspected_events", 5, 0, 5);
  //Run over all files in file1, presumably /hdfs/store/user/<etc>/ (must end with a /)                                                                       
  TString path = file1;
  TSystemDirectory sourceDir("hi",path);
  TList* fileList = sourceDir.GetListOfFiles();
  TIter next(fileList);
  TSystemFile* filename;
  int fileNumber = 0;
  int maxFiles = -1;
  while ((filename = (TSystemFile*)next()) && fileNumber >  maxFiles)
    {
      //if(fileNumber > 1)
	{
	  TString dataset = "Ntuple_";
	  TString  FullPathInputFile = (path+filename->GetName());
	  TString name = filename->GetName();

	  if(name.Contains(dataset))
	    {
	      //TFile *f_ins=new TFile(FullPathInputFile);
	      //TH1F *h_ins=(TH1F*) f_ins->Get("ins_Events");
	      //std::cout<<"FullPathInputFile:"<<FullPathInputFile<<std::endl;
	      chain->Add(FullPathInputFile);
	      //inspected_events->Add(h_ins);
	      //h_ins->Delete();
	      //f_ins->Close();
	    }	  
	}
	fileNumber++;
    }
  std::cout<<"********** for MC *********"<<std::endl;
  std::cout<<"All files added."<<std::endl;
  std::cout<<"Initializing chain."<<std::endl;
  Init(chain, isMC);
  BookHistos(file2);
 
}


skimm_et_2017::~skimm_et_2017()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   fileName->Write();
   fileName->Close();
}

Int_t skimm_et_2017::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skimm_et_2017::LoadTree(Long64_t entry)
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

void skimm_et_2017::Init(TChain *tree, string _isMC_)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
  TString isMC = TString(_isMC_);
  cout<<"from Init "<<isMC<<endl;
   // Set object pointer
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoUnCalibE = 0;
   phoUnCalibESigma = 0;
   phoCalibE = 0;
   phoCalibESigma = 0;
   phoCalibEt = 0;
   phoEnergyScale = 0;
   phoEnergySigma = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9Full5x5 = 0;
   phoHoverE = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIPhiFull5x5 = 0;
   phoSigmaIPhiIPhiFull5x5 = 0;
   phoPFChIso = 0;
   phoPFChWorstIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoIDMVA = 0;
   phoIDbit = 0;
   phoSeedTime = 0;
   phoSeedEnergy = 0;
   phoFiredSingleTrgs = 0;
   phoFiredDoubleTrgs = 0;
   phoFiredTripleTrgs = 0;
   phoFiredL1Trgs = 0;
   phoScale_up = 0;
   phoScale_dn = 0;
   phoScale_stat_up = 0;
   phoScale_stat_dn = 0;
   phoScale_syst_up = 0;
   phoScale_syst_dn = 0;
   phoScale_gain_up = 0;
   phoScale_gain_dn = 0;
   phoResol_up = 0;
   phoResol_dn = 0;
   phoResol_rho_up = 0;
   phoResol_rho_dn = 0;
   phoResol_phi_up = 0;
   phoResol_phi_dn = 0;
   jetPt = 0;
   jetE = 0;
   jetEta = 0;
   jetPhi = 0;
   jetRawPt = 0;
   jetRawE = 0;
   jetMt = 0;
   jetArea = 0;
   jetMass = 0;
   jetMaxDistance = 0;
   jetPhiPhiMoment = 0;
   jetConstituentEtaPhiSpread = 0;
   jetConstituentPtDistribution = 0;
   jetPileup = 0;
   jetID = 0;
   jetPUID = 0;
   jetPUFullID = 0;
   jetPartonID = 0;
   jetHadFlvr = 0;
   jetJECUnc = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetPhotonEnF = 0;
   jetElectronEnF = 0;
   jetMuonEnF = 0;
   jetChargedMuonEnF = 0;
   jetHFHAE = 0;
   jetHFEME = 0;
   jetNConst = 0;
   jetNConstituents = 0;
   jetNCharged = 0;
   jetNNeutral = 0;
   jetNChargedHad = 0;
   jetNNeutralHad = 0;
   jetNPhoton = 0;
   jetNElectron = 0;
   jetNMuon = 0;
   jetCSV2BJetTags = 0;
   jetDeepCSVTags_b = 0;
   jetDeepCSVTags_bb = 0;
   jetDeepCSVTags_c = 0;
   jetDeepCSVTags_udsg = 0;
   jetDeepFlavour_b = 0;
   jetDeepFlavour_bb = 0;
   jetDeepFlavour_lepb = 0;
   jetDeepFlavour_c = 0;
   jetDeepFlavour_uds = 0;
   jetDeepFlavour_g = 0;
   jetetaWidth = 0;
   jetphiWidth = 0;
   jetConstPt = 0;
   jetConstEt = 0;
   jetConstEta = 0;
   jetConstPhi = 0;
   jetConstPdgId = 0;
   jetGenJetE = 0;
   jetGenJetPt = 0;
   jetGenJetEta = 0;
   jetGenJetPhi = 0;
   jetGenPartonID = 0;
   jetGenE = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenPhi = 0;
   jetGenPartonMomID = 0;
   jetP4Smear = 0;
   jetP4SmearUp = 0;
   jetP4SmearDo = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9Full5x5 = 0;
   eleE = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleD0 = 0;
   eleDz = 0;
   eleSIP = 0;
   eleUnCalibE = 0;
   eleUnCalibESigma = 0;
   eleCalibEecalonly = 0;
   eleCalibE = 0;
   eleCalibESigma = 0;
   eleCalibEt = 0;
   eleCalibEtSigma = 0;
   eleEnergyScale = 0;
   eleEnergySigma = 0;
   eleSCRawE = 0;
   eleSCE = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   eleEoverP = 0;
   eleEoverPout = 0;
   eleEoverPInv = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eledEtaAtCalo = 0;
   eledEtaseedAtVtx = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   elePFClusEcalIso = 0;
   elePFClusHcalIso = 0;
   eleFiredSingleTrgs = 0;
   eleFiredDoubleTrgs = 0;
   eleFiredL1Trgs = 0;
   eleHEEPID = 0;
   eleMVAIsoID = 0;
   eleMVAnoIsoID = 0;
   eleIDbit = 0;
   eleTrkdxy = 0;
   eleKFHits = 0;
   eleKFChi2 = 0;
   eleGSFChi2 = 0;
   eleScale_up = 0;
   eleScale_dn = 0;
   eleScale_stat_up = 0;
   eleScale_stat_dn = 0;
   eleScale_syst_up = 0;
   eleScale_syst_dn = 0;
   eleScale_gain_up = 0;
   eleScale_gain_dn = 0;
   eleResol_up = 0;
   eleResol_dn = 0;
   eleResol_rho_up = 0;
   eleResol_rho_dn = 0;
   eleResol_phi_up = 0;
   eleResol_phi_dn = 0;
   muPt = 0;
   muE = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIDbit = 0;
   muD0 = 0;
   muDz = 0;
   muSIP = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muStations = 0;
   muMatches = 0;
   muTrkQuality = 0;
   muInnervalidFraction = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   muBestTrkType = 0;
   muFiredTrgs = 0;
   muFiredL1Trgs = 0;
   tau_Pt = 0;
   tau_Et = 0;
   tau_Eta = 0;
   tau_Phi = 0;
   tau_Charge = 0;
   tau_DecayMode = 0;
   tau_decayModeFindingNewDMs = 0;
   tau_P = 0;
   tau_Vz = 0;
   tau_Energy = 0;
   tau_Mass = 0;
   tau_Dxy = 0;
   tau_ZImpact = 0;
   tau_byCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   tau_chargedIsoPtSum = 0;
   tau_neutralIsoPtSum = 0;
   tau_neutralIsoPtSumWeight = 0;
   tau_footprintCorrection = 0;
   tau_photonPtSumOutsideSignalCone = 0;
   tau_puCorrPtSum = 0;
   tau_NumSignalPFChargedHadrCands = 0;
   tau_NumSignalPFNeutrHadrCands = 0;
   tau_NumSignalPFGammaCands = 0;
   tau_NumSignalPFCands = 0;
   tau_NumIsolationPFChargedHadrCands = 0;
   tau_NumIsolationPFNeutrHadrCands = 0;
   tau_NumIsolationPFGammaCands = 0;
   tau_NumIsolationPFCands = 0;
   tau_LeadChargedHadronEta = 0;
   tau_LeadChargedHadronPhi = 0;
   tau_LeadChargedHadronPt = 0;
   tau_LeadChargedHadron_dz = 0;
   tau_LeadChargedHadron_dxy = 0;
   tau_IDbits = 0;
   tau_byIsolationMVArun2017v2DBoldDMwLTraw2017 = 0;
   tau_byVVVLooseDeepTau2017v2p1VSjet = 0;
   tau_byVVLooseDeepTau2017v2p1VSjet = 0;
   tau_byVLooseDeepTau2017v2p1VSjet = 0;
   tau_byLooseDeepTau2017v2p1VSjet = 0;
   tau_byMediumDeepTau2017v2p1VSjet = 0;
   tau_byTightDeepTau2017v2p1VSjet = 0;
   tau_byVTightDeepTau2017v2p1VSjet = 0;
   tau_byVVTightDeepTau2017v2p1VSjet = 0;
   tau_byVVVLooseDeepTau2017v2p1VSe = 0;
   tau_byVVLooseDeepTau2017v2p1VSe = 0;
   tau_byVLooseDeepTau2017v2p1VSe = 0;
   tau_byLooseDeepTau2017v2p1VSe = 0;
   tau_byMediumDeepTau2017v2p1VSe = 0;
   tau_byTightDeepTau2017v2p1VSe = 0;
   tau_byVTightDeepTau2017v2p1VSe = 0;
   tau_byVVTightDeepTau2017v2p1VSe = 0;
   tau_byVLooseDeepTau2017v2p1VSmu = 0;
   tau_byLooseDeepTau2017v2p1VSmu = 0;
   tau_byMediumDeepTau2017v2p1VSmu = 0;
   tau_byTightDeepTau2017v2p1VSmu = 0;
   tauFiredTrgs = 0;
   tauFiredL1Trgs = 0;
   pdf = 0;
   pdfSystWeight = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   mcPID = 0;
   mcVtx = 0;
   mcVty = 0;
   mcVtz = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcStatus = 0;
   mcStatusFlag = 0;
   mcIndex = 0;
   mcDaughterPID = 0;
   mcCharge = 0;
   mcMotherPID = 0;
   mcMotherIndex = 0;
   mcMotherStatus = 0;
   mcDaughterStatus = 0;
   mcDaughterList = 0;
   mcTauDecayMode = 0;
   genMatch2 = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX);
   fChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY);
   fChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ);
   fChain->SetBranchAddress("vtxNtrks", &vtxNtrks, &b_vtxNtrks);
   fChain->SetBranchAddress("vtx_isFake", &vtx_isFake, &b_vtx_isFake);
   fChain->SetBranchAddress("vtx_ndof", &vtx_ndof, &b_vtx_ndof);
   fChain->SetBranchAddress("vtx_rho", &vtx_rho, &b_vtx_rho);
   fChain->SetBranchAddress("isGoodVtx", &isGoodVtx, &b_isGoodVtx);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
   fChain->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
   fChain->SetBranchAddress("HLTEleMuXRejectedByPS", &HLTEleMuXRejectedByPS, &b_HLTEleMuXRejectedByPS);
   fChain->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   fChain->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   fChain->SetBranchAddress("HLTPhoRejectedByPS", &HLTPhoRejectedByPS, &b_HLTPhoRejectedByPS);
   fChain->SetBranchAddress("HLTTau", &HLTTau, &b_HLTTau);
   fChain->SetBranchAddress("HLTTauIsPrescaled", &HLTTauIsPrescaled, &b_HLTTauIsPrescaled);
   fChain->SetBranchAddress("HLTTauRejectedByPS", &HLTTauRejectedByPS, &b_HLTTauRejectedByPS);
   fChain->SetBranchAddress("HLTMet", &HLTMet, &b_HLTMet);
   fChain->SetBranchAddress("HLTMetIsPrescaled", &HLTMetIsPrescaled, &b_HLTMetIsPrescaled);
   fChain->SetBranchAddress("HLTMetRejectedByPS", &HLTMetRejectedByPS, &b_HLTMetRejectedByPS);
   fChain->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
   fChain->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
   fChain->SetBranchAddress("HLTJetRejectedByPS", &HLTJetRejectedByPS, &b_HLTJetRejectedByPS);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoUnCalibE", &phoUnCalibE, &b_phoUnCalibE);
   fChain->SetBranchAddress("phoUnCalibESigma", &phoUnCalibESigma, &b_phoUnCalibESigma);
   fChain->SetBranchAddress("phoCalibE", &phoCalibE, &b_phoCalibE);
   fChain->SetBranchAddress("phoCalibESigma", &phoCalibESigma, &b_phoCalibESigma);
   fChain->SetBranchAddress("phoCalibEt", &phoCalibEt, &b_phoCalibEt);
   fChain->SetBranchAddress("phoEnergyScale", &phoEnergyScale, &b_phoEnergyScale);
   fChain->SetBranchAddress("phoEnergySigma", &phoEnergySigma, &b_phoEnergySigma);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedEnergy", &phoSeedEnergy, &b_phoSeedEnergy);
   fChain->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   fChain->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   fChain->SetBranchAddress("phoFiredTripleTrgs", &phoFiredTripleTrgs, &b_phoFiredTripleTrgs);
   fChain->SetBranchAddress("phoFiredL1Trgs", &phoFiredL1Trgs, &b_phoFiredL1Trgs);
   fChain->SetBranchAddress("phoScale_up", &phoScale_up, &b_phoScale_up);
   fChain->SetBranchAddress("phoScale_dn", &phoScale_dn, &b_phoScale_dn);
   fChain->SetBranchAddress("phoScale_stat_up", &phoScale_stat_up, &b_phoScale_stat_up);
   fChain->SetBranchAddress("phoScale_stat_dn", &phoScale_stat_dn, &b_phoScale_stat_dn);
   fChain->SetBranchAddress("phoScale_syst_up", &phoScale_syst_up, &b_phoScale_syst_up);
   fChain->SetBranchAddress("phoScale_syst_dn", &phoScale_syst_dn, &b_phoScale_syst_dn);
   fChain->SetBranchAddress("phoScale_gain_up", &phoScale_gain_up, &b_phoScale_gain_up);
   fChain->SetBranchAddress("phoScale_gain_dn", &phoScale_gain_dn, &b_phoScale_gain_dn);
   fChain->SetBranchAddress("phoResol_up", &phoResol_up, &b_phoResol_up);
   fChain->SetBranchAddress("phoResol_dn", &phoResol_dn, &b_phoResol_dn);
   fChain->SetBranchAddress("phoResol_rho_up", &phoResol_rho_up, &b_phoResol_rho_up);
   fChain->SetBranchAddress("phoResol_rho_dn", &phoResol_rho_dn, &b_phoResol_rho_dn);
   fChain->SetBranchAddress("phoResol_phi_up", &phoResol_phi_up, &b_phoResol_phi_up);
   fChain->SetBranchAddress("phoResol_phi_dn", &phoResol_phi_dn, &b_phoResol_phi_dn);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetE", &jetE, &b_jetE);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawE", &jetRawE, &b_jetRawE);
   fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetMaxDistance", &jetMaxDistance, &b_jetMaxDistance);
   fChain->SetBranchAddress("jetPhiPhiMoment", &jetPhiPhiMoment, &b_jetPhiPhiMoment);
   fChain->SetBranchAddress("jetConstituentEtaPhiSpread", &jetConstituentEtaPhiSpread, &b_jetConstituentEtaPhiSpread);
   fChain->SetBranchAddress("jetConstituentPtDistribution", &jetConstituentPtDistribution, &b_jetConstituentPtDistribution);
   fChain->SetBranchAddress("jetPileup", &jetPileup, &b_jetPileup);
   fChain->SetBranchAddress("jetID", &jetID, &b_jetID);
   fChain->SetBranchAddress("jetPUID", &jetPUID, &b_jetPUID);
   fChain->SetBranchAddress("jetPUFullID", &jetPUFullID, &b_jetPUFullID);
   fChain->SetBranchAddress("jetPartonID", &jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetHadFlvr", &jetHadFlvr, &b_jetHadFlvr);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetPhotonEnF", &jetPhotonEnF, &b_jetPhotonEnF);
   fChain->SetBranchAddress("jetElectronEnF", &jetElectronEnF, &b_jetElectronEnF);
   fChain->SetBranchAddress("jetMuonEnF", &jetMuonEnF, &b_jetMuonEnF);
   fChain->SetBranchAddress("jetChargedMuonEnF", &jetChargedMuonEnF, &b_jetChargedMuonEnF);
   fChain->SetBranchAddress("jetHFHAE", &jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", &jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConst", &jetNConst, &b_jetNConst);
   fChain->SetBranchAddress("jetNConstituents", &jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetNCharged", &jetNCharged, &b_jetNCharged);
   fChain->SetBranchAddress("jetNNeutral", &jetNNeutral, &b_jetNNeutral);
   fChain->SetBranchAddress("jetNChargedHad", &jetNChargedHad, &b_jetNChargedHad);
   fChain->SetBranchAddress("jetNNeutralHad", &jetNNeutralHad, &b_jetNNeutralHad);
   fChain->SetBranchAddress("jetNPhoton", &jetNPhoton, &b_jetNPhoton);
   fChain->SetBranchAddress("jetNElectron", &jetNElectron, &b_jetNElectron);
   fChain->SetBranchAddress("jetNMuon", &jetNMuon, &b_jetNMuon);
   fChain->SetBranchAddress("jetCSV2BJetTags", &jetCSV2BJetTags, &b_jetCSV2BJetTags);
   fChain->SetBranchAddress("jetDeepCSVTags_b", &jetDeepCSVTags_b, &b_jetDeepCSVTags_b);
   fChain->SetBranchAddress("jetDeepCSVTags_bb", &jetDeepCSVTags_bb, &b_jetDeepCSVTags_bb);
   fChain->SetBranchAddress("jetDeepCSVTags_c", &jetDeepCSVTags_c, &b_jetDeepCSVTags_c);
   fChain->SetBranchAddress("jetDeepCSVTags_udsg", &jetDeepCSVTags_udsg, &b_jetDeepCSVTags_udsg);
   fChain->SetBranchAddress("jetDeepFlavour_b", &jetDeepFlavour_b, &b_jetDeepFlavour_b);
   fChain->SetBranchAddress("jetDeepFlavour_bb", &jetDeepFlavour_bb, &b_jetDeepFlavour_bb);
   fChain->SetBranchAddress("jetDeepFlavour_lepb", &jetDeepFlavour_lepb, &b_jetDeepFlavour_lepb);
   fChain->SetBranchAddress("jetDeepFlavour_c", &jetDeepFlavour_c, &b_jetDeepFlavour_c);
   fChain->SetBranchAddress("jetDeepFlavour_uds", &jetDeepFlavour_uds, &b_jetDeepFlavour_uds);
   fChain->SetBranchAddress("jetDeepFlavour_g", &jetDeepFlavour_g, &b_jetDeepFlavour_g);
   fChain->SetBranchAddress("jetetaWidth", &jetetaWidth, &b_jetetaWidth);
   fChain->SetBranchAddress("jetphiWidth", &jetphiWidth, &b_jetphiWidth);
   fChain->SetBranchAddress("jetConstPt", &jetConstPt, &b_jetConstPt);
   fChain->SetBranchAddress("jetConstEt", &jetConstEt, &b_jetConstEt);
   fChain->SetBranchAddress("jetConstEta", &jetConstEta, &b_jetConstEta);
   fChain->SetBranchAddress("jetConstPhi", &jetConstPhi, &b_jetConstPhi);
   fChain->SetBranchAddress("jetConstPdgId", &jetConstPdgId, &b_jetConstPdgId);
   
   if(isMC=="MC"){
   
   fChain->SetBranchAddress("jetGenJetE", &jetGenJetE, &b_jetGenJetE);
   fChain->SetBranchAddress("jetGenJetPt", &jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", &jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", &jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", &jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenE", &jetGenE, &b_jetGenE);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("jetGenPartonMomID", &jetGenPartonMomID, &b_jetGenPartonMomID);
   fChain->SetBranchAddress("jetP4Smear", &jetP4Smear, &b_jetP4Smear);
   fChain->SetBranchAddress("jetP4SmearUp", &jetP4SmearUp, &b_jetP4SmearUp);
   fChain->SetBranchAddress("jetP4SmearDo", &jetP4SmearDo, &b_jetP4SmearDo);
   }
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleE", &eleE, &b_eleE);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleSIP", &eleSIP, &b_eleSIP);
   fChain->SetBranchAddress("eleUnCalibE", &eleUnCalibE, &b_eleUnCalibE);
   fChain->SetBranchAddress("eleUnCalibESigma", &eleUnCalibESigma, &b_eleUnCalibESigma);
   fChain->SetBranchAddress("eleCalibEecalonly", &eleCalibEecalonly, &b_eleCalibEecalonly);
   fChain->SetBranchAddress("eleCalibE", &eleCalibE, &b_eleCalibE);
   fChain->SetBranchAddress("eleCalibESigma", &eleCalibESigma, &b_eleCalibESigma);
   fChain->SetBranchAddress("eleCalibEt", &eleCalibEt, &b_eleCalibEt);
   fChain->SetBranchAddress("eleCalibEtSigma", &eleCalibEtSigma, &b_eleCalibEtSigma);
   fChain->SetBranchAddress("eleEnergyScale", &eleEnergyScale, &b_eleEnergyScale);
   fChain->SetBranchAddress("eleEnergySigma", &eleEnergySigma, &b_eleEnergySigma);
   fChain->SetBranchAddress("eleSCRawE", &eleSCRawE, &b_eleSCRawE);
   fChain->SetBranchAddress("eleSCE", &eleSCE, &b_eleSCE);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eledEtaAtCalo", &eledEtaAtCalo, &b_eledEtaAtCalo);
   fChain->SetBranchAddress("eledEtaseedAtVtx", &eledEtaseedAtVtx, &b_eledEtaseedAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso, &b_elePFClusEcalIso);
   fChain->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso, &b_elePFClusHcalIso);
   fChain->SetBranchAddress("eleFiredSingleTrgs", &eleFiredSingleTrgs, &b_eleFiredSingleTrgs);
   fChain->SetBranchAddress("eleFiredDoubleTrgs", &eleFiredDoubleTrgs, &b_eleFiredDoubleTrgs);
   fChain->SetBranchAddress("eleFiredL1Trgs", &eleFiredL1Trgs, &b_eleFiredL1Trgs);
   fChain->SetBranchAddress("eleHEEPID", &eleHEEPID, &b_eleHEEPID);
   fChain->SetBranchAddress("eleMVAIsoID", &eleMVAIsoID, &b_eleMVAIsoID);
   fChain->SetBranchAddress("eleMVAnoIsoID", &eleMVAnoIsoID, &b_eleMVAnoIsoID);
   fChain->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
   fChain->SetBranchAddress("eleTrkdxy", &eleTrkdxy, &b_eleTrkdxy);
   fChain->SetBranchAddress("eleKFHits", &eleKFHits, &b_eleKFHits);
   fChain->SetBranchAddress("eleKFChi2", &eleKFChi2, &b_eleKFChi2);
   fChain->SetBranchAddress("eleGSFChi2", &eleGSFChi2, &b_eleGSFChi2);
   fChain->SetBranchAddress("eleScale_up", &eleScale_up, &b_eleScale_up);
   fChain->SetBranchAddress("eleScale_dn", &eleScale_dn, &b_eleScale_dn);
   fChain->SetBranchAddress("eleScale_stat_up", &eleScale_stat_up, &b_eleScale_stat_up);
   fChain->SetBranchAddress("eleScale_stat_dn", &eleScale_stat_dn, &b_eleScale_stat_dn);
   fChain->SetBranchAddress("eleScale_syst_up", &eleScale_syst_up, &b_eleScale_syst_up);
   fChain->SetBranchAddress("eleScale_syst_dn", &eleScale_syst_dn, &b_eleScale_syst_dn);
   fChain->SetBranchAddress("eleScale_gain_up", &eleScale_gain_up, &b_eleScale_gain_up);
   fChain->SetBranchAddress("eleScale_gain_dn", &eleScale_gain_dn, &b_eleScale_gain_dn);
   fChain->SetBranchAddress("eleResol_up", &eleResol_up, &b_eleResol_up);
   fChain->SetBranchAddress("eleResol_dn", &eleResol_dn, &b_eleResol_dn);
   fChain->SetBranchAddress("eleResol_rho_up", &eleResol_rho_up, &b_eleResol_rho_up);
   fChain->SetBranchAddress("eleResol_rho_dn", &eleResol_rho_dn, &b_eleResol_rho_dn);
   fChain->SetBranchAddress("eleResol_phi_up", &eleResol_phi_up, &b_eleResol_phi_up);
   fChain->SetBranchAddress("eleResol_phi_dn", &eleResol_phi_dn, &b_eleResol_phi_dn);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muE", &muE, &b_muE);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIDbit", &muIDbit, &b_muIDbit);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muSIP", &muSIP, &b_muSIP);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("muBestTrkType", &muBestTrkType, &b_muBestTrkType);
   fChain->SetBranchAddress("muFiredTrgs", &muFiredTrgs, &b_muFiredTrgs);
   fChain->SetBranchAddress("muFiredL1Trgs", &muFiredL1Trgs, &b_muFiredL1Trgs);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("tau_Pt", &tau_Pt, &b_tau_Pt);
   fChain->SetBranchAddress("tau_Et", &tau_Et, &b_tau_Et);
   fChain->SetBranchAddress("tau_Eta", &tau_Eta, &b_tau_Eta);
   fChain->SetBranchAddress("tau_Phi", &tau_Phi, &b_tau_Phi);
   fChain->SetBranchAddress("tau_Charge", &tau_Charge, &b_tau_Charge);
   fChain->SetBranchAddress("tau_DecayMode", &tau_DecayMode, &b_tau_DecayMode);
   fChain->SetBranchAddress("tau_decayModeFindingNewDMs", &tau_decayModeFindingNewDMs, &b_tau_decayModeFindingNewDMs);
   fChain->SetBranchAddress("tau_P", &tau_P, &b_tau_P);
   fChain->SetBranchAddress("tau_Vz", &tau_Vz, &b_tau_Vz);
   fChain->SetBranchAddress("tau_Energy", &tau_Energy, &b_tau_Energy);
   fChain->SetBranchAddress("tau_Mass", &tau_Mass, &b_tau_Mass);
   fChain->SetBranchAddress("tau_Dxy", &tau_Dxy, &b_tau_Dxy);
   fChain->SetBranchAddress("tau_ZImpact", &tau_ZImpact, &b_tau_ZImpact);
   fChain->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("tau_chargedIsoPtSum", &tau_chargedIsoPtSum, &b_tau_chargedIsoPtSum);
   fChain->SetBranchAddress("tau_neutralIsoPtSum", &tau_neutralIsoPtSum, &b_tau_neutralIsoPtSum);
   fChain->SetBranchAddress("tau_neutralIsoPtSumWeight", &tau_neutralIsoPtSumWeight, &b_tau_neutralIsoPtSumWeight);
   fChain->SetBranchAddress("tau_footprintCorrection", &tau_footprintCorrection, &b_tau_footprintCorrection);
   fChain->SetBranchAddress("tau_photonPtSumOutsideSignalCone", &tau_photonPtSumOutsideSignalCone, &b_tau_photonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("tau_puCorrPtSum", &tau_puCorrPtSum, &b_tau_puCorrPtSum);
   fChain->SetBranchAddress("tau_NumSignalPFChargedHadrCands", &tau_NumSignalPFChargedHadrCands, &b_tau_NumSignalPFChargedHadrCands);
   fChain->SetBranchAddress("tau_NumSignalPFNeutrHadrCands", &tau_NumSignalPFNeutrHadrCands, &b_tau_NumSignalPFNeutrHadrCands);
   fChain->SetBranchAddress("tau_NumSignalPFGammaCands", &tau_NumSignalPFGammaCands, &b_tau_NumSignalPFGammaCands);
   fChain->SetBranchAddress("tau_NumSignalPFCands", &tau_NumSignalPFCands, &b_tau_NumSignalPFCands);
   fChain->SetBranchAddress("tau_NumIsolationPFChargedHadrCands", &tau_NumIsolationPFChargedHadrCands, &b_tau_NumIsolationPFChargedHadrCands);
   fChain->SetBranchAddress("tau_NumIsolationPFNeutrHadrCands", &tau_NumIsolationPFNeutrHadrCands, &b_tau_NumIsolationPFNeutrHadrCands);
   fChain->SetBranchAddress("tau_NumIsolationPFGammaCands", &tau_NumIsolationPFGammaCands, &b_tau_NumIsolationPFGammaCands);
   fChain->SetBranchAddress("tau_NumIsolationPFCands", &tau_NumIsolationPFCands, &b_tau_NumIsolationPFCands);
   fChain->SetBranchAddress("tau_LeadChargedHadronEta", &tau_LeadChargedHadronEta, &b_tau_LeadChargedHadronEta);
   fChain->SetBranchAddress("tau_LeadChargedHadronPhi", &tau_LeadChargedHadronPhi, &b_tau_LeadChargedHadronPhi);
   fChain->SetBranchAddress("tau_LeadChargedHadronPt", &tau_LeadChargedHadronPt, &b_tau_LeadChargedHadronPt);
   fChain->SetBranchAddress("tau_LeadChargedHadron_dz", &tau_LeadChargedHadron_dz, &b_tau_LeadChargedHadron_dz);
   fChain->SetBranchAddress("tau_LeadChargedHadron_dxy", &tau_LeadChargedHadron_dxy, &b_tau_LeadChargedHadron_dxy);
   fChain->SetBranchAddress("tau_IDbits", &tau_IDbits, &b_tau_IDbits);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMwLTraw2017", &tau_byIsolationMVArun2017v2DBoldDMwLTraw2017, &b_tau_byIsolationMVArun2017v2DBoldDMwLTraw2017);
   fChain->SetBranchAddress("tau_byVVVLooseDeepTau2017v2p1VSjet", &tau_byVVVLooseDeepTau2017v2p1VSjet, &b_tau_byVVVLooseDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byVVLooseDeepTau2017v2p1VSjet", &tau_byVVLooseDeepTau2017v2p1VSjet, &b_tau_byVVLooseDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byVLooseDeepTau2017v2p1VSjet", &tau_byVLooseDeepTau2017v2p1VSjet, &b_tau_byVLooseDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byLooseDeepTau2017v2p1VSjet", &tau_byLooseDeepTau2017v2p1VSjet, &b_tau_byLooseDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byMediumDeepTau2017v2p1VSjet", &tau_byMediumDeepTau2017v2p1VSjet, &b_tau_byMediumDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byTightDeepTau2017v2p1VSjet", &tau_byTightDeepTau2017v2p1VSjet, &b_tau_byTightDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byVTightDeepTau2017v2p1VSjet", &tau_byVTightDeepTau2017v2p1VSjet, &b_tau_byVTightDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byVVTightDeepTau2017v2p1VSjet", &tau_byVVTightDeepTau2017v2p1VSjet, &b_tau_byVVTightDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byVVVLooseDeepTau2017v2p1VSe", &tau_byVVVLooseDeepTau2017v2p1VSe, &b_tau_byVVVLooseDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVVLooseDeepTau2017v2p1VSe", &tau_byVVLooseDeepTau2017v2p1VSe, &b_tau_byVVLooseDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVLooseDeepTau2017v2p1VSe", &tau_byVLooseDeepTau2017v2p1VSe, &b_tau_byVLooseDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byLooseDeepTau2017v2p1VSe", &tau_byLooseDeepTau2017v2p1VSe, &b_tau_byLooseDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byMediumDeepTau2017v2p1VSe", &tau_byMediumDeepTau2017v2p1VSe, &b_tau_byMediumDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byTightDeepTau2017v2p1VSe", &tau_byTightDeepTau2017v2p1VSe, &b_tau_byTightDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVTightDeepTau2017v2p1VSe", &tau_byVTightDeepTau2017v2p1VSe, &b_tau_byVTightDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVVTightDeepTau2017v2p1VSe", &tau_byVVTightDeepTau2017v2p1VSe, &b_tau_byVVTightDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byVLooseDeepTau2017v2p1VSmu", &tau_byVLooseDeepTau2017v2p1VSmu, &b_tau_byVLooseDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("tau_byLooseDeepTau2017v2p1VSmu", &tau_byLooseDeepTau2017v2p1VSmu, &b_tau_byLooseDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("tau_byMediumDeepTau2017v2p1VSmu", &tau_byMediumDeepTau2017v2p1VSmu, &b_tau_byMediumDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("tau_byTightDeepTau2017v2p1VSmu", &tau_byTightDeepTau2017v2p1VSmu, &b_tau_byTightDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("tauFiredTrgs", &tauFiredTrgs, &b_tauFiredTrgs);
   fChain->SetBranchAddress("tauFiredL1Trgs", &tauFiredL1Trgs, &b_tauFiredL1Trgs);

   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("caloMET", &caloMET, &b_caloMET);
   fChain->SetBranchAddress("caloMETPhi", &caloMETPhi, &b_caloMETPhi);
   fChain->SetBranchAddress("caloMETsumEt", &caloMETsumEt, &b_caloMETsumEt);
   fChain->SetBranchAddress("pfMETCorr", &pfMETCorr, &b_pfMETCorr);
   fChain->SetBranchAddress("pfMETPhiCorr", &pfMETPhiCorr, &b_pfMETPhiCorr);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
   fChain->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
   fChain->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
   fChain->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
   fChain->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
   fChain->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
   fChain->SetBranchAddress("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp, &b_pfMETPhi_T1JESUp);
   fChain->SetBranchAddress("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo, &b_pfMETPhi_T1JESDo);
   fChain->SetBranchAddress("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp, &b_pfMETPhi_T1UESUp);
   fChain->SetBranchAddress("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo, &b_pfMETPhi_T1UESDo);

   if(isMC=="MC"){
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   
   fChain->SetBranchAddress("pdf", &pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("pdfWeight", &pdfWeight, &b_pdfWeight);
   fChain->SetBranchAddress("pdfSystWeight", &pdfSystWeight, &b_pdfSystWeight);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", &mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcVty", &mcVty, &b_mcVty);
   fChain->SetBranchAddress("mcVtz", &mcVtz, &b_mcVtz);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag);
   fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("mcDaughterPID", &mcDaughterPID, &b_mcDaughterPID);
   fChain->SetBranchAddress("mcCharge", &mcCharge, &b_mcCharge);
   fChain->SetBranchAddress("mcMotherPID", &mcMotherPID, &b_mcMotherPID);
   fChain->SetBranchAddress("mcMotherIndex", &mcMotherIndex, &b_mcMotherIndex);
   fChain->SetBranchAddress("mcMotherStatus", &mcMotherStatus, &b_mcMotherStatus);
   fChain->SetBranchAddress("mcDaughterStatus", &mcDaughterStatus, &b_mcDaughterStatus);
   fChain->SetBranchAddress("mcDaughterList", &mcDaughterList, &b_mcDaughterList);
   fChain->SetBranchAddress("mcTauDecayMode", &mcTauDecayMode, &b_mcTauDecayMode);
   fChain->SetBranchAddress("genMatch2", &genMatch2, &b_genMatch2); 
   }

   Notify();
}

Bool_t skimm_et_2017::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skimm_et_2017::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skimm_et_2017::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef skimm_et_2017_cxx
