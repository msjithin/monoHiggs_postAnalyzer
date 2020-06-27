
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
//#include "myevent.h"
//#include "LinkDef.h"
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"

using namespace std;


double test_sf(double muPt , double muEta, int tauDM, bool fakeBkg , bool isMC, bool debug)
{
  double rv_sf=1.0;
  double sf_muID = 1.0;
  double recoMuonPt=0.0;
  if (muPt < 120)
    recoMuonPt=muPt;
  else
    recoMuonPt = 119;

  TFile *f_muIDSF= TFile::Open("sf_files/RunBCDEF_SF_ID.root", "READ");
  TH2F *h_muIDSF=(TH2F*) f_muIDSF->Get("NUM_MediumID_DEN_genTracks_pt_abseta");
  
  sf_muID = h_muIDSF->GetBinContent(h_muIDSF->GetXaxis()->FindBin(recoMuonPt),h_muIDSF->GetYaxis()->FindBin(abs(muEta)));
  
  //if(debug)cout<<" setting up other files ..."<<endl;
  TFile *f_pileup = new TFile("sf_files/RootFiles/pileup/PU_Central.root");
  TH1F* h_pileup = (TH1F*)f_pileup->Get("pileup");
  
  TFile *f_muIDSF= TFile::Open("sf_files/RunBCDEF_SF_ID.root", "READ");
  //TFile *f_muIDSF= new TFile("sf_files/RootFiles/muon/2017_RunBCDEF_SF_ID.root");
  TH2F *h_muIDSF=(TH2F*) f_muIDSF->Get("NUM_MediumID_DEN_genTracks_pt_abseta");
  
  TFile *f_muIsoSF= TFile::Open("sf_files/RunBCDEF_SF_ISO.root", "READ");
  //TFile *f_muIsoSF= new TFile("sf_files/RootFiles/muon/2017_RunBCDEF_SF_ISO.root");
  TH2F *h_muIsoSF=(TH2F*) f_muIsoSF->Get("NUM_TightRelIso_DEN_MediumID_pt_abseta");

  TFile *f_muTrgSF= new TFile("sf_files/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root");
  TH2F *h_muTrgSF=(TH2F*) f_muTrgSF->Get("IsoMu27_PtEtaBins/pt_abseta_ratio");
  
  TFile *f_tauidSF = new TFile("sf_files/TauIDSFs/data/TauID_SF_dm_DeepTau2017v2p1VSjet_2017ReReco.root");
  TH1F *h_tauidSF_m = (TH1F*)f_tauidSF->Get("Medium");
  TH1F *h_tauidSF_vvvl = (TH1F*)f_tauidSF->Get("VVVLoose");
  
  /* TFile *f_tauidSF = new TFile("sf_files/TauIDSFs/data/TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco.root"); */
  /* TF1 *fn_tauIDSF_m = (TF1*) f_tauidSF->Get("Medium_cent"); */
  /* TF1 *fn_tauIDSF_vvl = (TF1*) f_tauidSF->Get("VVLoose_cent"); */

  TFile *f_tauesSF = new TFile("sf_files/TauIDSFs/data/TauES_dm_DeepTau2017v2p1VSjet_2017ReReco.root");
  TH1F *h_tauesSF = (TH1F*)f_tauesSF->Get("tes");
  
  TFile *f_tauFakeMuSF = new TFile("sf_files/TauIDSFs/data/TauID_SF_eta_DeepTau2017v2p1VSmu_2017ReReco.root");
  TH1F *h_tauFakeMuSF = (TH1F*)f_tauFakeMuSF->Get("Tight");
  
  TFile *f_tauFakeEleSF = new TFile("sf_files/TauIDSFs/data/TauID_SF_eta_DeepTau2017v2p1VSe_2017ReReco.root");
  TH1F *h_tauFakeEleSF = (TH1F*)f_tauFakeEleSF->Get("VVLoose");
  
  TFile *f_taufesSF = TFile::Open("sf_files/TauIDSFs/data/TauFES_eta-dm_DeepTau2017v2p1VSe_2017ReReco.root");
  TGraph *h_taufesSF = (TGraph*) f_taufesSF->Get("fes");
  //if(debug)cout<<" sf files opened successfully..."<<endl;

  rv_sf = sf_muID;
  return rv_sf;
  
}

#endif	/* _JETVETO_H */

