//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 14 01:21:53 2018 by ROOT version 6.10/09
// from TTree EventTree/Event data (tag V08_00_26_03)
// found on file: /hdfs/store/user/jmadhusu/MonoHiggs_MC2017/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50/180603_140815/0000/ggtree_mc_1.root
//////////////////////////////////////////////////////////

#ifndef mutau_analyzer_h
#define mutau_analyzer_h

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
#include "makeHisto.h"
// Header file for the classes stored in the TTree if any.
#include "vector"
#include "TString.h"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
//#include "makeHisto.h"

using namespace std;

class mutau_analyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   //TTree *outTree=new TTree("eventTree","eventTree");
   TFile *fileName;
   TTree *tree = new TTree("myTree", "A ROOT tree");
   TH1F* h_nEvents ;
   Int_t reco_mu_index=-1; Int_t reco_tau_index=-1;
   int nBackupTriggerEvents, nBTMediumEvents, nBTMediumHLTsinglePhoEvents, nEffPhoptden, nEffPhoptnum, nEffMETden, nEffMETnum;
   double nHLTPassed,n_eventWeight, nSingleTrgPassed, nGoodMuonPassed, nElectronPtPassed, nGoodTauPassed, nTauPtPassed,numberOfEvents,nMETPassed, nDPhiPassed, nqcdden,nqcdnum, nMETFiltersPassed,nLeptonVetoPassed,nPassedBjetVeto,nNoisyCrystals,nDPhiJetMETPassed, nGoodMuTauPassed,nDeltaRPassed,nPassedThirdLepVeto,nMtPassed,  nPassedHiggsPtcut, nPassedVisibleMasscut, nPassedMETcut, nFinal_afterSelections, nGoodMuonPassed_qcd, nGoodTauPassed_qcd, nGoodMuTauPassed_qcd, nDeltaRPassed_qcd, nPassedThirdLepVeto_qcd,nPassedBjetVeto_qcd, nPassedHiggsPtcut_qcd, nPassedVisibleMasscut_qcd, nPassedMETcut_qcd, nFinal_afterSelections_qcd ,nPasstottrmass, nPassJetsSelection,nL1PrefiringPassed, nMuonDzPassed, nMuonD0Passed, nMuonIdPassed, nMuonIsoPassed, nTauIsoPassed, nTauDecayModePassed, nTauRejectionPassed;
  double nPassedSkimmed, nEtaCutsPassed, nMetfiltersPassed, nHTTTriger, nAntileptondiscriminators, nLeptonseparation, nbtaggingVetos, nTauPtcut,nGenJetRemoval, nOppositeChargePassed,nMediumDeepTauID;

   
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

   //   mutau_analyzer(TTree *tree=0);
   mutau_analyzer(const char* file1, const char* file2, string isMC);
   virtual ~mutau_analyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree, string isMC);
   //virtual void     Init(TTree *tree);
   virtual void     Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void BookHistos(const char* file1, const char* file2);
   virtual void fillHistos(int histoNumber, double event_weight,int higgs_index);
   virtual double DeltaPhi(double phi1, double phi2);
   virtual vector<int> getMuCand(double pt, double eta);
   virtual vector<int> getTauCand(double pt, double eta);
   virtual vector<int> found_higgs();
   virtual vector<int> found_muon();
   virtual vector<int> found_electron();
   virtual vector<int> found_tau();
   virtual vector<int> found_tauh();
   virtual vector<int> found_tauNeu();
   virtual bool skimming_Htt();
   virtual vector<int> skimmed_Mu(); 
   virtual vector<int> skimmed_Tau();
   virtual int gen_matching();
   virtual int thirdLeptonVeto();
   virtual double dR(int mu_index, int tau_index);
   virtual double delta_R(float phi1, float eta1, float phi2, float eta2);
   virtual float TMass_F(float LepPt, float LepPhi ,float met, float metPhi);
   virtual float TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met);   
   virtual float VisMass_F(TLorentzVector a, TLorentzVector b);
   virtual float pTvecsum_F(float pt1, float pt2, float phi1, float phi2);
   //   virtual bool electron_pass(int pho_index, float elePtCut);
   //virtual bool relIso(int ele_index);
   virtual bool passBjetVeto();
   virtual void fillHist( string histNumber, int muIndex, int tauIndex, float event_weight);
   virtual void fillHist( string histNumber, TLorentzVector muonP4, TLorentzVector tauP4, int muIndex, int tauIndex, float event_weight);
   virtual vector<int> getGenMuPlus();
   virtual vector<int> getGenMuMinus();


};

#endif

#ifdef mutau_analyzer_cxx
mutau_analyzer::mutau_analyzer(const char* file1, const char* file2, string isMC)
{
  TChain *chain = new TChain("eventTree");
  
  TString  FullPathInputFile = file1;
  chain->Add(FullPathInputFile);
  
  std::cout<<"All files added."<<std::endl;
  std::cout<<"Initializing chain."<<std::endl;
  Init(chain, isMC);
  BookHistos(file1, file2);

  //inspected_events->Write();
}


mutau_analyzer::~mutau_analyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   /*   map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
   map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
   for (; iMap1 != jMap1; ++iMap1)
     nplot1(iMap1->first)->Write();   
   */
   fileName->Write();
   /*   map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
   map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
   for (; iMap1 != jMap1; ++iMap1)
     nplot1(iMap1->first)->Write();
   */
   tree->Write();
   fileName->Close();
}

Int_t mutau_analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mutau_analyzer::LoadTree(Long64_t entry)
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

void mutau_analyzer::Init(TChain *tree, string _isMC_)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
    nBackupTriggerEvents = nBTMediumEvents = nBTMediumHLTsinglePhoEvents = nEffPhoptden = nEffPhoptnum = nEffMETden = nEffMETnum = 0;
    
    nHLTPassed = n_eventWeight = nSingleTrgPassed = nGoodMuonPassed = nElectronPtPassed = nGoodTauPassed = nTauPtPassed= numberOfEvents = nMETPassed = nDPhiPassed = nqcdden= nqcdnum=nMETFiltersPassed= nLeptonVetoPassed=nPassedBjetVeto=nNoisyCrystals=nDPhiJetMETPassed= nGoodMuTauPassed = nDeltaRPassed= nPassedThirdLepVeto=nMtPassed=nPassedHiggsPtcut=nPassedVisibleMasscut=nPassedMETcut=nFinal_afterSelections=nGoodMuonPassed_qcd=nGoodTauPassed_qcd=nGoodMuTauPassed_qcd=nDeltaRPassed_qcd=nPassedThirdLepVeto_qcd=nPassedBjetVeto_qcd=nPassedHiggsPtcut_qcd=nPassedVisibleMasscut_qcd=nPassedMETcut_qcd=nFinal_afterSelections_qcd=nPasstottrmass= nPassJetsSelection=nL1PrefiringPassed=nMuonDzPassed=nMuonD0Passed=nMuonIdPassed=nMuonIsoPassed=nTauIsoPassed=nTauDecayModePassed=nTauRejectionPassed=0;

    nPassedSkimmed=nEtaCutsPassed=nMetfiltersPassed=nHTTTriger=nAntileptondiscriminators=nLeptonseparation=nbtaggingVetos=nTauPtcut=nGenJetRemoval=nOppositeChargePassed=nMediumDeepTauID=0;
   

  TString isMC = TString(_isMC_);
  cout<<"from Init "<<isMC<<endl;
   // Set object pointer
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
   mcPID = 0;
   mcPt = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcStatus = 0;
   mcStatusFlag = 0;
   mcIndex = 0;
   mcDaughterPID = 0;
   mcCharge = 0;
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
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   if(isMC=="MC"){
     fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
     fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
     fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
     fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
     fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
     fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
     fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
     fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
     fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
     fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag);
     fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
     fChain->SetBranchAddress("mcDaughterPID", &mcDaughterPID, &b_mcDaughterPID);
     fChain->SetBranchAddress("mcCharge", &mcCharge, &b_mcCharge);
     fChain->SetBranchAddress("genMatch2", &genMatch2, &b_genMatch2);
   }
   Notify();
}

Bool_t mutau_analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mutau_analyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mutau_analyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mutau_analyzer_cxx
