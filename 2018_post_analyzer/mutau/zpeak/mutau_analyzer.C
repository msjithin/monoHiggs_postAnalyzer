/////mutau_analyzer.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom mutau_analyzer analyze
//
//To run, assuming this is compiled to an executable named 'analyze':
//$ ./analyze /hdfs/store/user/jmadhusu/LatestNtuples/ /afs/hep.wisc.edu/user/ms/CMSSW_9_4_4/src/2017_analysis/etau/output.root -1 10000
//./analyze /hdfs/store/user/jmadhusu/MonoHiggs_MC2017/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/crab_ZZZ/180603_185329/0000/ /afs/hep.wisc.edu/user/ms/CMSSW_9_4_4/src/2017_analysis/analyzer/output.root -1 10000
//Runs over every event in the folder LatestNtuples, reporting progress every 10000 events
//and storing the resulting histograms in the file output.root.
//
//To plot, for example, single photon trigger efficiency as a function of photon pt:
//$ root -l
//root[0] TFile *f = new TFile("output.root");
//root[1] TGraphAsymmErrors *efficiency = new TGraphAsymmErrors((TH1F*)f->Get("Photon_Et_300_2"),(TH1F*)f->Get("Photon_Et_300_1"));
//root[2] efficiency->Draw("AP")
//root[3] efficiency->SetTitle("Single photon trigger efficiency")
//root[4] efficiency->GetXaxis()->SetTitle("Photon p_{T}")
//root[5] efficiency->Draw("AP")
//
#define mutau_analyzer_cxx
#include "mutau_analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
#include <TLorentzVector.h>
#include "makeHisto.h" 
#include "/afs/hep.wisc.edu/home/ms/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/post_analysis/mutau/skimm_submit/skimmed/sf_files/roCorr_Run2_v3/RoccoR.cc"

using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{
  TStopwatch sw;
  sw.Start();
  
  myMap1 = new map<string, TH1F*>();
  //myMap2 = new map<string, TH2F*>();
  std::string SampleName = argv[7];
  std::string isMC  = argv[6];
  std::string outputfile  = argv[2];
  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
    {
      std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
      return 1;
    }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
    {
      std::cout<<"Please enter a valid value for reportEvery (parameter 4) "<<std::endl;
      return 1;
    }
  //std::string SampleName = argv[5];
  
  mutau_analyzer t(argv[1],argv[2], isMC);
  t.Loop(maxEvents,reportEvery, SampleName , isMC);
  //delete myMap1;
  cout<<" Outpt written to "<<outputfile<<endl;
  sw.Stop();
  sw.Print();
  return 0;
}

void mutau_analyzer::Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
{
  
  int nTotal;
  nTotal = 0;
  int report_=0;
  int report_test=0;
  int nInspected;
  nInspected = 0;
  double nInspected_genWeighted;  
  nInspected_genWeighted = 0.0; 
  bool debug=false;  
  double netWeight = 1.0;
  double afterSF1=0;
  double afterSF2=0;     
  double afterSF3=0;     
  double afterSF4=0;     

  if (fChain == 0) return;
  int genMatching; 
  int thirdLeptonIndex=-1;
  std::vector<int> muonGen;
  std::vector<int> muCand;        muCand.clear();
  std::vector<int> tauCand;       tauCand.clear();
  
  std::vector<int> higgsCand;     higgsCand.clear();
  std::vector<int> reco_mu;       reco_mu.clear(); 
  std::vector<int> reco_tau;      reco_tau.clear(); 
  
  std::vector<int> eleGenCand;    eleGenCand.clear();
  std::vector<int> muGenCand;     muGenCand.clear();   
  std::vector<int> tauGenCand;    tauGenCand.clear();
  std::vector<int> tauhGenCand;   tauhGenCand.clear();
  std::vector<int> tauNeuGenCand; tauNeuGenCand.clear();
  
  std::vector<int> skimmedMuon;   skimmedMuon.clear(); 
  std::vector<int> skimmedTau;    skimmedTau.clear();
  
  TString sample = TString(SampleName);
  
  int nHiggs = 0;
   int nHToMuTau = 0;
   int found_mt = 0;
   int muCand_1=0; int muCand_2=0;int muCand_3=0;
   int tauCand_1=0; int tauCand_2=0;int tauCand_3=0;

   bool fill_hist = false;
   bool isMC = false;
   if( _isMC_=="MC" ) { isMC=true; fill_hist=true; }
   else if ( _isMC_=="DATA" ) { isMC=false; fill_hist=false; }
  
   Double_t  Pt_Bins[26]={0.0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
   Double_t  Pt_Bins_highPt[21]={100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};

   TH1F* h_cutflow=new TH1F("cutflow", "cutflow", 10, 0, 10); h_cutflow->Sumw2();
   TH1F* h_cutflow_n=new TH1F("cutflow_n", "cutflow_n", 8, 0, 8);h_cutflow_n->Sumw2();
   
   TH1F* h_cutflow_Htt=new TH1F("cutflow_Htt", "cutflow_Htt", 11, 0, 11); h_cutflow_Htt->Sumw2();

   TLorentzVector myMomTau, myTauh,  myNeu, myHiggs; 
   TFile *f_kfactors = TFile::Open("sf_files/kfactors.root", "READ");
   TH1D *ewkCorrection = (TH1D*)f_kfactors->Get("EWKcorr/W");
   TH1D *NNLOCorrection = (TH1D*)f_kfactors->Get("WJets_LO/inv_pt");
   bool found_Wjet_sample=false;
   //W1JetsToLNu, 
   
   
   if ( sample.Contains("WJetsToLNu") ||
	sample.Contains("W1JetsToLNu") ||
	sample.Contains("W2JetsToLNu") ||
	sample.Contains("W3JetsToLNu") ||
	sample.Contains("W4JetsToLNu") 	) {
     found_Wjet_sample=true;
     cout<<"******************wjet sample found"<<endl;
   }
   
   RoccoR  rc("sf_files/roCorr_Run2_v3/RoccoR2018.txt"); 
   TFile* f_muon_id = TFile::Open("sf_files/muon/2018_RunABCD_SF_ID.root", "READ");
   TFile* f_muon_iso= TFile::Open("sf_files/muon/2018_RunABCD_SF_ISO.root", "READ");
   TFile* f_ele_reco = TFile::Open("sf_files/electron/2018_egammaEffi_txt_EGM2D_updatedAll.root", "READ");
   TFile* f_ele_id_noiso = TFile::Open("sf_files/electron/2018_ElectronWPVeto_Fall17V2.root", "READ");
   TFile* f_ele_id_tight = TFile::Open("sf_files/electron/2018_ElectronTight.root", "READ");
   TFile* f_trigger_sf_1 = TFile::Open("sf_files/trigger/EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root", "READ");
   TFile* f_trigger_sf_2 = TFile::Open("sf_files/trigger/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root", "READ");
   
   TH2F* h_muon_id  = (TH2F*)((TH2F*)f_muon_id->Get("NUM_MediumID_DEN_TrackerMuons_pt_abseta"))->Clone(TString("NUM_MediumID_DEN_TrackerMuons_pt_abseta"));
   TH2F* h_muon_iso  = (TH2F*)((TH2F*)f_muon_iso->Get("NUM_TightRelIso_DEN_MediumID_pt_abseta"))->Clone(TString("NUM_TightRelIso_DEN_MediumID_pt_abseta"));
   TH2F* h_muon_trg_1  = (TH2F*)((TH2F*)f_trigger_sf_1->Get("IsoMu24_PtEtaBins/pt_abseta_ratio"))->Clone(TString("pt_abseta_ratio"));
   TH2F* h_muon_trg_2  = (TH2F*)((TH2F*)f_trigger_sf_2->Get("IsoMu24_PtEtaBins/pt_abseta_ratio"))->Clone(TString("pt_abseta_ratio"));
   TH2F* h_ele_id_noiso  = (TH2F*)((TH2F*)f_ele_id_noiso->Get("EGamma_SF2D"))->Clone(TString("EGamma_SF2D"));
   h_muon_id->SetDirectory(0);
   h_muon_iso->SetDirectory(0);
   h_muon_trg_1->SetDirectory(0);
   h_muon_trg_2->SetDirectory(0);
   h_ele_id_noiso->SetDirectory(0);


   f_muon_id->Close();
   f_muon_iso->Close();
   f_ele_reco->Close();
   f_ele_id_noiso->Close();
   f_ele_id_tight->Close();
   f_trigger_sf_1->Close();
   f_trigger_sf_2->Close();
   
   TLorentzVector muonP4;
   TLorentzVector muonPlusP4; TLorentzVector muonMinusP4;
   TLorentzVector tauP4;

   Long64_t nentries = fChain->GetEntries();
   if ( isMC==true ) std::cout<<".... MC file ..... "<<std::endl;
   else  std::cout<<".... DATA file ..... "<<std::endl;

   std::cout<<"Coming in: "<<std::endl;
   std::cout<<"nentries:"<<nentries<<std::endl;
   //Look at up to maxEvents events, or all if maxEvents == -1.
   Long64_t nentriesToCheck = nentries;
   if (maxEvents != -1LL && nentries > maxEvents)
     nentriesToCheck = maxEvents;
   nTotal = nentriesToCheck;
   Long64_t nbytes = 0, nb = 0;
   
   std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
   //TStopwatch sw;
   //sw.Start();
   
   for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
     {
       
       //event_.clear();
       //event_info.clear();
       higgsCand.clear();
       muCand.clear();
       tauCand.clear();
       reco_mu.clear(); 
       reco_tau.clear();  
       eleGenCand.clear();   
       muGenCand.clear();
       tauGenCand.clear();
       tauhGenCand.clear();
       tauNeuGenCand.clear();
       skimmedMuon.clear();
       skimmedTau.clear();
       muonGen.clear();
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       double inspected_event_weight = 1.0; 
       if(isMC)	 fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight/fabs(genWeight) : inspected_event_weight = 0.0;
       nInspected_genWeighted += inspected_event_weight;  
       nInspected += 1; 
       //h_insEvents->SetBinContent(1, nInspected_genWeighted);
       //=1.0 for real data
       double event_weight=1.0;
       int report_i=0;
       numberOfEvents+=event_weight;
       event_weight=inspected_event_weight;
       //cout<<jentry<< "  nMC : "<<nMC<<endl;
       //cout<<jentry<< "  genMatching : "<<genMatching<<endl;
       double muRC_sf = 1.0; 
       double muon_id_sf = 1.0;
       double muon_iso_sf = 1.0;
       double ele_reco_sf = 1.0;
       double ele_id_sf = 1.0;
       double trigger_sf = 1.0;
       
       double EWK_corrected_weight=1.0;
       double NNLO_weight = 1.0;
       double kfactor = 1.0;
       
       int bosonPID;
       double bosonPt;
       bool Wfound = false;
       if (isMC)
	 {       //check which mc particle is W boson                                                                                               
	   for(int i=0; i<nMC;i++){
	     if((*mcPID)[i] == 24){
	       Wfound=true;
	       bosonPID = (*mcPID)[i];
	       bosonPt = (*mcPt)[i];
	     }
	   }
	 }
	 
       if (isMC) genMatching = gen_matching();
       else genMatching = 0;
       //if(debug)cout<<"this worked Line 202"<<endl;
       //// adding skimming
       
       if(isMC) event_weight=inspected_event_weight;
       else event_weight=1.0;
       int leading_muon = -1; float leading_mPt=0;
       int leading_tau = -1;  float leading_tPt=0;
       ////// reco selection begin
       if(debug)cout<<"this worked Line 278"<<endl;
       if(metFilters==0)
	 {
	   //fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed+=event_weight;
	   //if(fill_hist==true )//fillHistos(0,event_weight, higgsCand[0]);
	   if( HLTEleMuX>>21&1 == 1 ||
	       HLTEleMuX>>60&1 == 1    )
	     {
	       nSingleTrgPassed+=event_weight;
	       //if(fill_hist==true )//fillHistos(1,event_weight, higgsCand[0]);
	       if(debug)cout<<"this worked Line 289"<<endl;

	       muCand = getMuCand(20,2.4);  ///// muons selected 
	       if( muCand.size() >0 ) 
		 { 
		   nGoodMuonPassed+=event_weight;
		   //cout<<"this worked Line 403"<<endl;
		   
		   if (found_Wjet_sample==true )
		     {
		       EWK_corrected_weight = 1.0*(ewkCorrection->GetBinContent(ewkCorrection->GetXaxis()->FindBin(bosonPt)));
		       NNLO_weight = 1.0*(NNLOCorrection->GetBinContent(NNLOCorrection->GetXaxis()->FindBin(bosonPt)));
		       if(EWK_corrected_weight!=0 && NNLO_weight!=0)
			 kfactor = (EWK_corrected_weight/NNLO_weight);
		       else
			 kfactor=1.21;
		       event_weight*=kfactor;
		     }
		   
		   std::vector<int> iMuPlus; iMuPlus.clear(); 
		   std::vector<int> iMuMinus; iMuMinus.clear();
		   for(int i=0; i<muCand.size(); i++){
		     if(muCharge->at(muCand[i]) < 0) iMuMinus.push_back(muCand[i]);
		     if(muCharge->at(muCand[i]) > 0) iMuPlus.push_back(muCand[i]);
		   }
		   if(iMuPlus.size()>0 && iMuMinus.size()>0)
		     {
		       if(debug)cout<<"this worked Line 316"<<endl;
		       muonPlusP4.SetPtEtaPhiE(muPt->at(iMuPlus[0]), muEta->at(iMuPlus[0]), muPhi->at(iMuPlus[0]), muE->at(iMuPlus[0]));
		       muonMinusP4.SetPtEtaPhiE(muPt->at(iMuMinus[0]), muEta->at(iMuMinus[0]), muPhi->at(iMuMinus[0]), muE->at(iMuMinus[0]));
		       
		       if(debug)cout<<"this worked Line 320"<<endl;
		       nGoodMuTauPassed+=event_weight;
		       double randomN = gRandom->Rndm();
		       
		       if(isMC==false) 
			 muRC_sf = rc.kScaleDT(muCharge->at(iMuPlus[0]), muPt->at(iMuPlus[0]), muEta->at(iMuPlus[0]), muPhi->at(iMuPlus[0]), 0, 0); //data
		       else {
			 if (genMatching==2 || genMatching==4){
			   muonGen=getGenMuPlus();
			   double mcGenPt=0;
			   if(muonGen.size()>0 )  mcGenPt= mcPt->at(muonGen[0]);
			   muRC_sf = rc.kSpreadMC(muCharge->at(iMuPlus[0]), muPt->at(iMuPlus[0]), muEta->at(iMuPlus[0]), muPhi->at(iMuPlus[0]), mcGenPt, 0, 0); 
			 }
			 else{
			   muRC_sf = rc.kSmearMC(muCharge->at(iMuPlus[0]), muPt->at(iMuPlus[0]), muEta->at(iMuPlus[0]), muPhi->at(iMuPlus[0]),muTrkLayers->at(iMuPlus[0]), randomN, 0, 0);
			 }
		       }
		       muonPlusP4 = muonPlusP4*muRC_sf;
		       afterSF1+=event_weight;
		       if(isMC==false)
                         muRC_sf = rc.kScaleDT(muCharge->at(iMuMinus[0]), muPt->at(iMuMinus[0]), muEta->at(iMuMinus[0]), muPhi->at(iMuMinus[0]), 0, 0); //data
                       else {
                         if (genMatching==2 || genMatching==4){
                           muonGen=getGenMuMinus();
			   double mcGenPt=0;
			   if(muonGen.size()>0 )  mcGenPt= mcPt->at(muonGen[0]);
                           muRC_sf = rc.kSpreadMC(muCharge->at(iMuMinus[0]), muPt->at(iMuMinus[0]), muEta->at(iMuMinus[0]), muPhi->at(iMuMinus[0]), mcGenPt, 0, 0);
                         }
                         else{
                           muRC_sf = rc.kSmearMC(muCharge->at(iMuMinus[0]), muPt->at(iMuMinus[0]), muEta->at(iMuMinus[0]), muPhi->at(iMuMinus[0]),muTrkLayers->at(iMuMinus[0]), randomN, 0, 0);
                         }
                       }
                       muonMinusP4 = muonMinusP4*muRC_sf;
                       afterSF1+=event_weight;
		       ////////////////////////////////
		       if(muonPlusP4.Pt()<120){
			 muon_id_sf  = h_muon_id->GetBinContent(h_muon_id->GetXaxis()->FindBin(muonPlusP4.Pt()), 
								h_muon_id->GetYaxis()->FindBin(fabs(muonPlusP4.Eta())));
			 muon_iso_sf = h_muon_iso->GetBinContent(h_muon_id->GetXaxis()->FindBin(muonPlusP4.Pt()),
								 h_muon_id->GetYaxis()->FindBin(fabs(muonPlusP4.Eta())));
		       }
		       else{
			 muon_id_sf  = h_muon_id->GetBinContent(h_muon_id->GetXaxis()->FindBin(119),
								h_muon_id->GetYaxis()->FindBin(fabs(muonPlusP4.Eta())));
			 muon_iso_sf = h_muon_iso->GetBinContent(h_muon_id->GetXaxis()->FindBin(119),
								 h_muon_id->GetYaxis()->FindBin(fabs(muonPlusP4.Eta())));
		       }
		       if (run < 316361)
			 trigger_sf  = h_muon_trg_1->GetBinContent(h_muon_id->GetXaxis()->FindBin(muonPlusP4.Pt()),
								   h_muon_id->GetYaxis()->FindBin(fabs(muonPlusP4.Eta())));
		       else if (run >= 316361)
			 trigger_sf  = h_muon_trg_2->GetBinContent(h_muon_id->GetXaxis()->FindBin(muonPlusP4.Pt()),
								   h_muon_id->GetYaxis()->FindBin(fabs(muonPlusP4.Eta())));
		       
		       if(isMC)event_weight = event_weight*muon_id_sf* muon_iso_sf* trigger_sf;
		       ////////////////////////////////
		       muon_id_sf=muon_iso_sf=trigger_sf=1.0;
		       if(muonMinusP4.Pt()<120){
                         muon_id_sf  = h_muon_id->GetBinContent(h_muon_id->GetXaxis()->FindBin(muonMinusP4.Pt()),
                                                                h_muon_id->GetYaxis()->FindBin(fabs(muonMinusP4.Eta())));
                         muon_iso_sf = h_muon_iso->GetBinContent(h_muon_id->GetXaxis()->FindBin(muonMinusP4.Pt()),
                                                                 h_muon_id->GetYaxis()->FindBin(fabs(muonMinusP4.Eta())));
                       }
                       else{
                         muon_id_sf  = h_muon_id->GetBinContent(h_muon_id->GetXaxis()->FindBin(119),
                                                                h_muon_id->GetYaxis()->FindBin(fabs(muonMinusP4.Eta())));
                         muon_iso_sf = h_muon_iso->GetBinContent(h_muon_id->GetXaxis()->FindBin(119),
                                                                 h_muon_id->GetYaxis()->FindBin(fabs(muonMinusP4.Eta())));
                       }
                       if (run < 316361)
                         trigger_sf  = h_muon_trg_1->GetBinContent(h_muon_id->GetXaxis()->FindBin(muonMinusP4.Pt()),
                                                                   h_muon_id->GetYaxis()->FindBin(fabs(muonMinusP4.Eta())));
                       else if (run >= 316361)
                         trigger_sf  = h_muon_trg_2->GetBinContent(h_muon_id->GetXaxis()->FindBin(muonMinusP4.Pt()),
                                                                   h_muon_id->GetYaxis()->FindBin(fabs(muonMinusP4.Eta())));

                       if(isMC)event_weight = event_weight*muon_id_sf* muon_iso_sf* trigger_sf;

		       afterSF2+=event_weight;
		       
		       if(debug) cout<<"before plot fill"<<endl;
		       plotFill("muPt_Plus_5",  muonPlusP4.Pt() , 40 , 0 , 200,  event_weight);
		       plotFill("muEta_Plus_5", muonPlusP4.Eta(), 30, -2.4, 2.4,  event_weight);
		       plotFill("muPhi_Plus_5", muonPlusP4.Phi(), 30, -3.14, 3.14,  event_weight);
		       plotFill("muDz_Plus_5",  muDz->at(iMuPlus[0]), 40, -0.5, 0.5,  event_weight);
		       plotFill("muD0_Plus_5",  muD0->at(iMuPlus[0]), 48, -0.06, 0.06,  event_weight);
		       plotFill("muonID_Plus_5", muIDbit->at(iMuPlus[0])>>1&1, 4, -2, 2,  event_weight); // muonID
		       plotFill("relMuIso_Plus_5", muIDbit->at(iMuPlus[0])>>8&1, 4, -2, 2,  event_weight);
		       plotFill("muCharge_Plus_5", muCharge->at(iMuPlus[0]), 8, -2, 2 ,  event_weight);
		       plotFill("muPt_Minus_5",  muonMinusP4.Pt() , 40 , 0 , 200,  event_weight);
                       plotFill("muEta_Minus_5", muonMinusP4.Eta(), 30, -2.4, 2.4,  event_weight);
                       plotFill("muPhi_Minus_5", muonMinusP4.Phi(), 30, -3.14, 3.14,  event_weight);
                       plotFill("muDz_Minus_5",  muDz->at(iMuMinus[0]), 40, -0.5, 0.5,  event_weight);
                       plotFill("muD0_Minus_5",  muD0->at(iMuMinus[0]), 48, -0.06, 0.06,  event_weight);
                       plotFill("muonID_Minus_5", muIDbit->at(iMuMinus[0])>>1&1, 4, -2, 2,  event_weight); // muonID
                       plotFill("relMuIso_Minus_5", muIDbit->at(iMuMinus[0])>>8&1, 4, -2, 2,  event_weight);
                       plotFill("muCharge_Minus_5", muCharge->at(iMuMinus[0]), 8, -2, 2 ,  event_weight);
		       float muMuMass=VisMass_F(muonPlusP4, muonMinusP4);
		       plotFill("muMuMass_5",  muMuMass , 50 , 0 , 200,  event_weight);
		     }
		 }
	     }
	 }
     
       report_test = nentriesToCheck/20;
       while (report_test>10)
	 {
	   report_test=report_test/10;
	   report_i++;
	 }
       reportEvery = report_test*pow(10,report_i);
       if (jentry%reportEvery == 0)
	 {
	   std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
	 }
     }
   
   

   std::cout.setf( std::ios::fixed, std:: ios::floatfield );
   if((nentriesToCheck-1)%reportEvery != 0)
     std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
   // sw.Stop();
   std::cout<<"All events checked."<<std::endl;
   std::cout<<"*******************************************"<<std::endl;
   std::cout<<"******************Jithin's original*************************"<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"Initial entries "<<numberOfEvents<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"Passing smikking "<<nPassedSkimmed<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"Inspected genWeightd "<<nInspected_genWeighted<<std::setw(10) <<std::right << "   % change= "<<(numberOfEvents-nInspected_genWeighted)*100/numberOfEvents<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"METFiltersPassed "<<nMETFiltersPassed<<std::setw(10) <<std::right << "   % change= "<<(nInspected_genWeighted-nMETFiltersPassed)*100/nInspected_genWeighted<<std::endl;

   std::cout<<std::setw(20) <<std::right <<"SingleTrgPassed "<<nSingleTrgPassed<<std::setw(10) <<std::right << "   % change= "<<(nMETFiltersPassed-nSingleTrgPassed)*100/nMETFiltersPassed<<std::endl;

   std::cout<<std::setw(20) <<std::right <<"GoodMuonPassed "<<nGoodMuonPassed<<std::setw(10) <<std::right << "   % change= "<<(nSingleTrgPassed-nGoodMuonPassed)*100/nSingleTrgPassed<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"GoodTauPassed "<<nGoodTauPassed<<std::setw(10) <<std::right << "   % change= "<<(nGoodMuonPassed-nGoodTauPassed)*100/nGoodMuonPassed<<std::endl;
   //   std::cout<<std::setw(20) <<std::right <<"TauIsoPassed "<<nTauIsoPassed<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nTauIsoPassed)*100/nGoodTauPassed<<std::endl;
   //std::cout<<std::setw(20) <<std::right <<"TauDecayModePassed "<<nTauDecayModePassed<<std::setw(10) <<std::right << "   % change= "<<(nTauIsoPassed-nTauDecayModePassed)*100/nTauIsoPassed<<std::endl;

   std::cout<<std::setw(20) <<std::right <<"opp charge "<<nGoodMuTauPassed<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nGoodMuTauPassed)*100/nGoodTauPassed<<std::endl;

   std::cout<<std::setw(20) <<std::right <<"after sf 1 "<<afterSF1<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nGoodMuTauPassed)*100/nGoodTauPassed<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"after sf 2 "<<afterSF2<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nGoodMuTauPassed)*100/nGoodTauPassed<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"after sf 3 "<<afterSF3<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nGoodMuTauPassed)*100/nGoodTauPassed<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"after sf 4 "<<afterSF4<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nGoodMuTauPassed)*100/nGoodTauPassed<<std::endl;

   
   std::cout<<std::setw(20) <<std::right <<"PassedThirdLepVeto "<<nPassedThirdLepVeto<<std::setw(10) <<std::right << "   % change= "<<(nGoodMuTauPassed-nPassedThirdLepVeto)*100/nGoodMuTauPassed<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"PassedBjetVeto "<<nPassedBjetVeto<<std::setw(10) <<std::right << "   % change= "<<(nPassedThirdLepVeto-nPassedBjetVeto)*100/nPassedThirdLepVeto<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"DeltaRPassed "<<nDeltaRPassed<<std::setw(10) <<std::right << "   % change= "<<(nPassedBjetVeto-nDeltaRPassed)*100/nPassedBjetVeto<<std::endl;


   std::cout<<std::setw(20) <<std::right <<"Total change :"<<(numberOfEvents-nDeltaRPassed)*100/numberOfEvents<<std::endl;
   std::cout<<"*******************************************"<<std::endl;
   std::cout<<"*******************************************"<<std::endl;
   std::cout<<std::setw(20) <<std::right <<"Number of events inspected: " << nInspected <<std::endl;
   std::cout<<std::setw(20) <<std::right << "Number of events inspected (minus negative gen. weights): " << nInspected_genWeighted << std::endl; 
   

   
   h_cutflow->SetBinContent(1,nInspected_genWeighted );
   //h_cutflow->SetBinContent(1, nHToMuTau );
   h_cutflow->SetBinContent(2, nSingleTrgPassed);
   h_cutflow->SetBinContent(3, nGoodMuonPassed);
   h_cutflow->SetBinContent(4, nGoodTauPassed);
   h_cutflow->SetBinContent(5, nTauIsoPassed);
   h_cutflow->SetBinContent(6, nTauDecayModePassed);
   h_cutflow->SetBinContent(7, nGoodMuTauPassed);
   h_cutflow->SetBinContent(8, nPassedThirdLepVeto);
   h_cutflow->SetBinContent(9, nPassedBjetVeto);
   h_cutflow->SetBinContent(10, nDeltaRPassed);
   
   h_cutflow_n->SetBinContent(1,nInspected_genWeighted );
   h_cutflow_n->SetBinContent(2, nSingleTrgPassed);
   h_cutflow_n->SetBinContent(3, nGoodMuonPassed);
   h_cutflow_n->SetBinContent(4, nGoodTauPassed);
   h_cutflow_n->SetBinContent(5, nGoodMuTauPassed);
   h_cutflow_n->SetBinContent(6, nPassedThirdLepVeto);
   h_cutflow_n->SetBinContent(7, nPassedBjetVeto);
   h_cutflow_n->SetBinContent(8, nDeltaRPassed);
     
   h_cutflow_Htt->SetBinContent( 1,nPassedSkimmed );
   h_cutflow_Htt->SetBinContent( 2,nEtaCutsPassed );
   h_cutflow_Htt->SetBinContent( 3,nMetfiltersPassed );
   h_cutflow_Htt->SetBinContent( 4,nHTTTriger );
   h_cutflow_Htt->SetBinContent( 5,nAntileptondiscriminators );
   h_cutflow_Htt->SetBinContent( 6,nLeptonseparation );
   h_cutflow_Htt->SetBinContent( 7,nbtaggingVetos );
   h_cutflow_Htt->SetBinContent( 8,nTauPtcut );
   h_cutflow_Htt->SetBinContent( 9,nGenJetRemoval );
   h_cutflow_Htt->SetBinContent(10,nOppositeChargePassed );
   h_cutflow_Htt->SetBinContent(11,nMediumDeepTauID );
   
   fileName->cd();
   map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
   map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
   for (; iMap1 != jMap1; ++iMap1)
     nplot1(iMap1->first)->Write();
}

void mutau_analyzer::BookHistos(const char* file1, const char* file2)
{
  TFile* file_in= new TFile(file1, "READ");
  fileName = new TFile(file2, "RECREATE");
  
  //makeOutputTree(tree);
  fileName->cd();
  h_nEvents = (TH1F*)((TH1F*)file_in->Get("nEvents"))->Clone(TString("nEvents"));
  file_in->Close();
  Float_t Pt_Bins[36]={0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
  
  Float_t MetBins[15]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180., 200, 300., 400., 600.0,800.0};
  Float_t TrMassBins[24]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180., 200, 220, 240,260,280,300.,320,340,360,380, 400., 600.0,800.0, 1000.0};
  
  /*
  //Set up the histos to be filled with method fillHistos
  for(int i=0; i<21; i++)
    {
      char ptbins[100];
      sprintf(ptbins, "_%d", i);
      std::string histname(ptbins);
      Double_t  Pt_Bins[26]={0.0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
      //h_dR[i] = new TH1F(("h_dR"+histname).c_str(),"h_dR",20,0,2.0);h_dR[i]->Sumw2();
      //h_HiggsPt[i]= new TH1F(("HiggsPt"+histname).c_str(),"HiggsPt", 100, 0.0, 1000.0);h_HiggsPt[i]->Sumw2();
      h_HiggsPt[i]= new TH1F(("Higgs_pt"+histname).c_str(),("Higgs_pt"+histname).c_str(), 25, Pt_Bins);h_HiggsPt[i]->Sumw2();
      //h_VisibleMass[i]= new TH1F(("VisibleMass"+histname).c_str(),"VisibleMass",20, 0, 200);h_VisibleMass[i]->Sumw2();      
    }
*/
}

//Fill the sequential histos at a particular spot in the sequence


void mutau_analyzer::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{
  
  //  h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);

}




//---------------------------------------------------                                                                                                                                
// get a electron candiate based on pt eta and isolation                                                                                                                               
//----------------------------------------------------                                                                                                                               

std::vector<int> mutau_analyzer::getMuCand(double muPtCut, double muEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over muons                                                                     
   for(int iMu=0;iMu<nMu;iMu++)
     {
       bool kinematic = false;
       bool muonID = false;
       bool muonIso =  false;
       bool trigger = false;
       if( muPt->at(iMu) > muPtCut  && fabs(muEta->at(iMu))< muEtaCut  && fabs(muDz->at(iMu)) < 0.2 && fabs(muD0->at(iMu))<0.045 ) kinematic = true;
       if( muIDbit->at(iMu)>>1&1==1 ) muonID = true;//|| muIDbit->at(iMu)>>8&1==1 || muIDbit->at(iMu)>>17&1==1  ) muonID = true;
       //float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
       //if( relMuIso < 0.15 ) muonIso = true;
       if( muIDbit->at(iMu)>>8&1==1 ) muonIso = true; // PFMedium isolation
       if(  (HLTEleMuX>>21&1 == 1 && muPt->at(iMu)>28) 
           || (HLTEleMuX>>60&1 == 1 && muPt->at(iMu)>25) 
           || (   HLTTau>>0&1  == 1 && muPt->at(iMu)>21 && muPt->at(iMu)<25 && muEta->at(iMu)<2.4 ) 
           )trigger=true;
       
       if( kinematic==true  &&  muonID==true &&  muonIso==true && trigger==true){
	 tmpCand.push_back(iMu);
       }                                                                                      
     }                                                                                       
  return tmpCand;
  
}

std::vector<int> mutau_analyzer::getTauCand(double tauPtCut, double tauEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over taus      
  for(int iTau=0;iTau<nTau;iTau++)
    {
      bool kinematic = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool mutau_separation=false;
      bool newDecayModeFinding=false;
      bool tau_reject=false;
      if( tau_Pt->at(iTau) > tauPtCut 
	  && fabs( tau_Eta->at(iTau))< tauEtaCut 
	  //&& fabs(tau_Charge->at(iTau))==1   
	  //&& fabs(tau_ZImpact->at(iTau)) < 200 
	  //&& fabs(tau_Dxy->at(iTau)) < 0.2 
	  )kinematic = true;
      if( tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1 ) tauIsolation=true; 
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 ||tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byTightDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;
      if( tau_decayModeFindingNewDMs->at(iTau)==1 ) newDecayModeFinding=true;
      

      if( kinematic==true    
	  && decayModeCut==true   
	  && tauIsolation==true 
	  && tau_reject==true   
	  && newDecayModeFinding==true
	  )
	{
	  tmpCand.push_back(iTau);
    	}                                                           
    }                                                                                       
  return tmpCand;
  
}

int mutau_analyzer::thirdLeptonVeto(){
  std::vector<int> tmpCand;
  tmpCand.clear();
  int thirdLepIndex = -1;
  bool thirdLepVeto=true;
  for(int iEle=0; iEle < nEle;iEle++)
    {
      bool kinematic = false;
      if( (*elePt)[iEle] > 10.0  && fabs((*eleEta)[iEle])< 2.5 && (*eleD0)[iEle] < 0.045 && (*eleDz)[iEle] < 0.2 ) kinematic = true;
      bool electronId =false;
      if( eleIDbit->at(iEle)>>8&1==1) electronId =true;
      bool relative_iso = false;
      float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
      if( relEleIso < 0.3 ) relative_iso = true;
      if(electronId==true && kinematic==true && relative_iso==true){
	tmpCand.push_back(iEle);
      }                                                                                         }                                                                                                                    
  if(tmpCand.size() > 0){ thirdLepIndex = tmpCand[0]; thirdLepVeto=false;}
  return thirdLepIndex;
  
}
                                                                                    

double mutau_analyzer::dR(int mu_index, int tau_index)
{
  double deltaeta = abs(muEta->at(mu_index) - tau_Eta->at(tau_index));
  double muonPhi = muPhi->at(mu_index);
  double tauPhi = tau_Phi->at(tau_index);

  double deltaphi = DeltaPhi(muonPhi, tauPhi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}

double mutau_analyzer::delta_R(float phi1, float eta1, float phi2, float eta2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar   = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}



double mutau_analyzer::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
}

float mutau_analyzer::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
    return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

float mutau_analyzer::TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float totalTMass = (a + b+ met).M();
  return totalTMass;
}


float mutau_analyzer::VisMass_F(TLorentzVector a, TLorentzVector b){
  float visibleMass = (a + b).M();
  return visibleMass;
}

float mutau_analyzer::pTvecsum_F(float pt1, float pt2, float phi1, float phi2) {
  float pt_vecSum = sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2), 2) + pow(pt1*sin(phi1) + pt2*sin(phi2), 2));
  return pt_vecSum;
}

bool mutau_analyzer::passBjetVeto()
{
  int tmpCand = 0;
  bool veto = true;
  bool foundBjet = false;
  for(int iJets=0; iJets<nJet ; iJets++){
    if(jetCSV2BJetTags->at(iJets) > 0.8838) tmpCand++;
    // CSV B jet tag for selecting bJets is medium WP (jetCSV2BJetTags > 0.8838.)
  }
  if(tmpCand>0){ veto = false; foundBjet = true; } 
  return veto;
}
std::vector<int> mutau_analyzer::found_higgs(){ 
  std::vector<int> tmpCand;    tmpCand.clear();   
  bool found_H=false;
  for(int i=0; i<nMC;i++){
    if (fabs(mcPID->at(i))==25  && mcStatus->at(i)==62 )
      {
	tmpCand.push_back(i);
      }
  }
  return tmpCand; 
}
std::vector<int> mutau_analyzer::found_muon(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int i=0; i<nMC;i++){
    if (fabs(mcPID->at(i))==13  )tmpCand.push_back(i); ///&& mcStatus->at(i)==1 
  }
  return tmpCand;
}
std::vector<int> mutau_analyzer::found_electron(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int i=0; i<nMC;i++){
    if (fabs(mcPID->at(i))==11 )tmpCand.push_back(i); ///&& mcStatus->at(i)==1
  }
  return tmpCand;
}

std::vector<int> mutau_analyzer::found_tau(){
  std::vector<int> tmpCand;    tmpCand.clear(); 
  for(int i=0; i<nMC;i++){
    if ( fabs(mcPID->at(i)) ==15 )tmpCand.push_back(i);
  }
  return tmpCand;
}
std::vector<int> mutau_analyzer::found_tauh(){
  std::vector<int> tmpCand;    tmpCand.clear();
  bool found_T=false;
  bool found_d=false;
  for(int i=0; i<nMC;i++){
    if ( fabs(mcPID->at(i))==15  ) found_T=true;
    if ( mcTauDecayMode->at(i)>>2&1==1 || mcTauDecayMode->at(i)>>3&1==1 || mcTauDecayMode->at(i)>>4&1==1 || mcTauDecayMode->at(i)>>5&1==1
	 || mcTauDecayMode->at(i)>>6&1==1 || mcTauDecayMode->at(i)>>7&1==1 || mcTauDecayMode->at(i)>>8&1==1 
	 || mcTauDecayMode->at(i)>>9&1==1 || mcTauDecayMode->at(i)>>10&1==1 || mcTauDecayMode->at(i)>>11&1==1) found_d=true;
    if (found_d==true )tmpCand.push_back(i);
  }
  return tmpCand;
}

std::vector<int> mutau_analyzer::found_tauNeu(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int i=0; i<nMC;i++){
    if (fabs(mcPID->at(i))==16)tmpCand.push_back(i);
  }
  return tmpCand;
}
bool mutau_analyzer::skimming_Htt(){
  bool tmpCand=false;
  bool muFound=false; bool tauFound=false; bool drCutPassed=false;
  std::vector<int> tmpMuCand;    tmpMuCand.clear();
  std::vector<int> tmpTauCand;    tmpTauCand.clear();

  for(int iMu=0; iMu<nMu;iMu++){
    if(fabs(muDz->at(iMu)) < 0.2 && 
       fabs(muD0->at(iMu))<0.045 && 
       muPt->at(iMu) > 19.5      &&
       fabs(muEta->at(iMu))< 2.4 &&
       muIDbit->at(iMu)>>8&1==1  
       ) {
      muFound=true;
      tmpMuCand.push_back(iMu);
    }
  }
  
  for(int iTau=0; iTau<nTau;iTau++){
    if( fabs(tau_ZImpact->at(iTau)) < 200 && 
	tau_Pt->at(iTau) > 29.5           &&
	fabs( tau_Eta->at(iTau))< 2.3     &&
	(tau_DecayMode->at(iTau) !=5 && tau_DecayMode->at(iTau)!=6) &&
	(tau_IDbits->at(iTau)>>2&1==1 || tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1 ) &&
	(tau_IDbits->at(iTau)>>4&1==1 || tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 || tau_IDbits->at(iTau)>>19&1==1)&&
	(tau_IDbits->at(iTau)>>14&1==1 || tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 )
	) {
      tauFound = true;
      tmpTauCand.push_back(iTau);
    }
  }
  
  for(int iMu=0; iMu<tmpMuCand.size();iMu++){
    for(int iTau=0; iTau<tmpTauCand.size();iTau++){
      double deltaR = dR(tmpMuCand[iMu], tmpTauCand[iTau]); 
      if(deltaR > 0.5 )
	drCutPassed=true;
    }
  }
  if(muFound==true && tauFound == true && drCutPassed==true)
    tmpCand=true;
  return tmpCand;
}

std::vector<int> mutau_analyzer::skimmed_Mu(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int iMu=0; iMu<nMu;iMu++){
    if(fabs(muDz->at(iMu)) < 0.2 &&
       fabs(muD0->at(iMu))<0.045 &&
       muPt->at(iMu) > 19.5      &&
       fabs(muEta->at(iMu))< 2.4 &&
       muIDbit->at(iMu)>>8&1==1
       ) {
      tmpCand.push_back(iMu);
    }
  }
  return tmpCand; 
}

std::vector<int> mutau_analyzer::skimmed_Tau(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int iTau=0; iTau<nTau;iTau++){
    if( fabs(tau_ZImpact->at(iTau)) < 200 &&
	tau_Pt->at(iTau) > 29.5           &&
	fabs( tau_Eta->at(iTau))< 2.3     &&
	(tau_DecayMode->at(iTau) !=5 && tau_DecayMode->at(iTau)!=6) &&
	(tau_IDbits->at(iTau)>>2&1==1 || tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1 ) &&
	(tau_IDbits->at(iTau)>>4&1==1 || tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 || tau_IDbits->at(iTau)>>19&1==1)&&
	(tau_IDbits->at(iTau)>>14&1==1 || tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 )
	) {
      tmpCand.push_back(iTau);
    }
  }
  return tmpCand;
}
int mutau_analyzer::gen_matching(){
  int tmpCand=-1;
  std::vector<int> tmpGenMatch;
  tmpGenMatch.clear();

  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>1&1==1 ) tmpGenMatch.push_back(1);
    if( genMatch2->at(imc)>>2&1==1 ) tmpGenMatch.push_back(2);
    if( genMatch2->at(imc)>>3&1==1 ) tmpGenMatch.push_back(3);
    if( genMatch2->at(imc)>>4&1==1 ) tmpGenMatch.push_back(4);
    if( genMatch2->at(imc)>>5&1==1 ) tmpGenMatch.push_back(5);
    if( genMatch2->at(imc)>>6&1==1 ) tmpGenMatch.push_back(6);
  }
  if(tmpGenMatch.size() >0 )
    tmpCand=tmpGenMatch[0];
  return tmpCand; 
}
std::vector<int> mutau_analyzer::getGenMuPlus(){
  std::vector<int> tmpCand;
  tmpCand.clear();
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>2&1==1 || genMatch2->at(imc)>>4&1==1 ) 
      if(mcCharge->at(imc) >0 )
	tmpCand.push_back(imc);
  }
  return tmpCand; 
}
std::vector<int> mutau_analyzer::getGenMuMinus(){
  std::vector<int> tmpCand;
  tmpCand.clear();
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>2&1==1 || genMatch2->at(imc)>>4&1==1 )
      if(mcCharge->at(imc) <0 )
        tmpCand.push_back(imc);
  }
  return tmpCand;
}


void mutau_analyzer::fillHist( string histNumber , int muIndex, int tauIndex, float event_weight){
  string hNumber = histNumber;
  plotFill("muPt_"+hNumber,  muPt->at(muIndex) , 40 , 0 , 200,  event_weight);
  plotFill("muEta_"+hNumber, muEta->at(muIndex), 30, -2.4, 2.4,  event_weight);
  plotFill("muPhi_"+hNumber, muPhi->at(muIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("muDz_"+hNumber,  muDz->at(muIndex), 20, -0.5, 0.5,  event_weight);
  plotFill("muD0_"+hNumber,  muD0->at(muIndex), 24, -0.06, 0.06,  event_weight);
  plotFill("muonID_"+hNumber,muIDbit->at(muIndex)>>1&1, 4, -2, 2,  event_weight); // muonID
  plotFill("relMuIso_"+hNumber, muIDbit->at(muIndex)>>8&1, 4, -2, 2,  event_weight);
  plotFill("muCharge_"+hNumber, muCharge->at(muIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 40 , 0 , 200,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 30, -2.3, 2.3,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byVVVLooseDeepTau2017v2p1VSe->at(tauIndex)==1, 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byTightDeepTau2017v2p1VSmu->at(tauIndex)==1, 8, -2, 2 ,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}
void mutau_analyzer::fillHist( string histNumber , TLorentzVector muonP4, TLorentzVector tauP4, int muIndex, int tauIndex, float event_weight){
  string hNumber = histNumber;
  plotFill("muPt_"+hNumber,  muonP4.Pt() , 40 , 0 , 200,  event_weight);
  plotFill("muEta_"+hNumber, muonP4.Eta(), 30, -2.4, 2.4,  event_weight);
  plotFill("muPhi_"+hNumber, muonP4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("muDz_"+hNumber,  muDz->at(muIndex), 20, -0.5, 0.5,  event_weight);
  plotFill("muD0_"+hNumber,  muD0->at(muIndex), 24, -0.06, 0.06,  event_weight);
  plotFill("muonID_"+hNumber,muIDbit->at(muIndex)>>1&1, 4, -2, 2,  event_weight); // muonID
  plotFill("relMuIso_"+hNumber, muIDbit->at(muIndex)>>8&1, 4, -2, 2,  event_weight);
  plotFill("muCharge_"+hNumber, muCharge->at(muIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 40 , 0 , 200,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 30, -2.3, 2.3,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byVVVLooseDeepTau2017v2p1VSe->at(tauIndex)==1, 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byTightDeepTau2017v2p1VSmu->at(tauIndex)==1, 8, -2, 2 ,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}
