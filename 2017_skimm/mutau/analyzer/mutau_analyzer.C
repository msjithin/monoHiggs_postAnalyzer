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
#include "sf_files/roCorr_Run2_v3/RoccoR.cc"

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
  double numberOfEvents = 0;
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
  int genMatching=0; 
  int thirdLeptonIndex=-1;
  std::vector<int> muonGen;       muonGen.clear();
  std::vector<int> muCand;        muCand.clear();
  std::vector<int> mu2Cand;       mu2Cand.clear();
  std::vector<int> tauCand;       tauCand.clear();
  std::vector<int> jetCand;       jetCand.clear();
  std::vector<int> higgsCand;     higgsCand.clear();
  std::vector<int> reco_mu;       reco_mu.clear(); 
  std::vector<int> reco_mu2;      reco_mu2.clear();
  std::vector<int> reco_tau;      reco_tau.clear(); 
    
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

   //TH1F* h_cutflow=new TH1F("cutflow", "cutflow", 10, 0, 10); h_cutflow->Sumw2();
   TH1F* h_cutflow_n=new TH1F("cutflow_n", "cutflow_n", 8, 0, 8);h_cutflow_n->Sumw2();
   TH1F* h_cutflow_n_fr=new TH1F("cutflow_n_fr", "cutflow_n_fr", 8, 0, 8);h_cutflow_n_fr->Sumw2();
   TH1F* h_cutflow_n_dyll=new TH1F("cutflow_n_dyll", "cutflow_n_dyll", 8, 0, 8);h_cutflow_n_dyll->Sumw2();
   TH1F* h_cutflow_n_dyll_fr=new TH1F("cutflow_n_dyll_fr", "cutflow_n_dyll_fr", 8, 0, 8);h_cutflow_n_dyll_fr->Sumw2();
   //TH1F* h_cutflow_Htt=new TH1F("cutflow_Htt", "cutflow_Htt", 11, 0, 11); h_cutflow_Htt->Sumw2();

   TLorentzVector myMomTau, myTauh,  myNeu, myHiggs; 
   bool found_Wjet_sample=false;
   bool found_DYjet_sample=false;
   int PID =0;
   if ( sample.Contains("WJetsToLNu") ||
	sample.Contains("W1JetsToLNu") ||
	sample.Contains("W2JetsToLNu") ||
	sample.Contains("W3JetsToLNu") ||
	sample.Contains("W4JetsToLNu") 	) {
     found_Wjet_sample=true;
     PID = 24;
     cout<<"****************** wjet sample found"<<endl;
   }
   if ( sample.Contains("DYJetsToLL") ||
	sample.Contains("DY1JetsToLL") ||
	sample.Contains("DY2JetsToLL") ||
	sample.Contains("DY3JetsToLL") ||
	sample.Contains("DY4JetsToLL")  ) {
     found_DYjet_sample=true;
     PID = 23;
     cout<<"****************** dyjet sample found"<<endl;
   }
   if(debug)cout<<" setting up kFactor files ..."<<endl;
   TH1F *NLO_QCD_EWK,*NLO_EWK,*NLO_QCD,*NNLO_QCD;
   TFile* f_nnlo_qcd = TFile::Open("sf_files/RootFiles/theory/lindert_qcd_nnlo_sf.root");
   TFile* f_nlo_qcd  = TFile::Open("sf_files/RootFiles/theory/2017_gen_v_pt_qcd_sf.root");
   TFile* f_qcd_ewk;
   if ( found_Wjet_sample ) {
     f_qcd_ewk = TFile::Open("sf_files/RootFiles/theory/merged_kfactors_wjets.root");
     NLO_QCD_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_qcd_ewk");
     NLO_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_ewk");
     NLO_QCD = (TH1F*)f_nlo_qcd->Get("wjet_dress_monojet");
     NNLO_QCD = (TH1F*)f_nnlo_qcd->Get("evj");
     
   } else if ( found_DYjet_sample ) {
     f_qcd_ewk = TFile::Open("sf_files/RootFiles/theory/merged_kfactors_zjets.root");
     NLO_QCD_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_qcd_ewk");
     NLO_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_ewk");
     f_nlo_qcd = TFile::Open("sf_files/RootFiles/theory/kfac_dy_filter.root");
     NLO_QCD = (TH1F*)f_nlo_qcd->Get("kfac_dy_filter");
     NNLO_QCD = (TH1F*)f_nnlo_qcd->Get("eej");
     
   }
   
   // if(debug)cout<<" setting up other files ..."<<endl;
   RoccoR  rc("sf_files/roCorr_Run2_v3/RoccoR2017.txt"); 
   TLorentzVector muonP4;
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
       if(debug) cout<<"event "<<jentry<<endl;
       muCand.clear(); mu2Cand.clear();
       tauCand.clear();
       jetCand.clear();
       reco_mu.clear(); 
       reco_mu2.clear();
       reco_tau.clear();  
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
       double weight=1.0;
       double muRC_sf = 1.0; double randomN = gRandom->Rndm();
       double sf_tauID = 1.0; 

       double pileup_sf = 1.0;
       double kfactor = 1.0;
       double nlo_ewk = 1.0;
       double nlo_qcd_binned = 1.0;
       double nlo_qcd = 1.0;
       double nnlo_qcd =1.0;
       int bosonPID;
       double bosonPt=0.0;
       bool Wfound = false;
       bool passSingleTriggerPaths=false;
       bool passCrossTrigger=false;
       int report_i=0;
       bool muTriggerFilterMatch=false;
       bool tauTriggerFilterMatch=false;
       bool muTau_selector=false;
       bool muMu_selector=false;
       numberOfEvents+=weight;
       weight=inspected_event_weight;
       if (isMC) genMatching = gen_matching();
       else genMatching = 0;
       //cout<<"genMatch = "<<genMatching<<endl;
       if(debug)cout<<"this worked Line 314"<<endl;
       
       if(isMC) weight=inspected_event_weight;
       else weight=1.0;
       int leading_muon = -1; float leading_mPt=0;
       int leading_tau = -1;  float leading_tPt=0;
       
       if(isMC)
	 pileup_sf = h_pileup->GetBinContent(h_pileup->GetXaxis()->FindBin(puTrue->at(0)));
       weight = weight*pileup_sf;
       if( isGoodVtx==false ) continue;
       //if( noisyJet2017()==true ) continue;
       //if( found_DYjet_sample && !(genMatching<5))
       //	 continue;

       if( found_DYjet_sample && hasGenTau())
	 muTau_selector=true;
       else if( found_DYjet_sample && !hasGenTau() )
	 muTau_selector=false;
       else if ( !found_DYjet_sample )
	 muTau_selector=true;

       if( found_DYjet_sample && genMatching<5 )
	 muMu_selector=true;
       else
	 muMu_selector=false;

       /////Trigger bit selection
       if(HLTEleMuX>>21&1 == 1 || HLTEleMuX>>60&1 == 1 )
	 passSingleTriggerPaths=true;
       if( HLTTau>>0&1 == 1  )
	 passCrossTrigger=true;
       ////
       if(isMC && (found_Wjet_sample || found_DYjet_sample)){
	 if(debug)cout<<"check which mc particle is W boson"<<endl;
	 for(int i=0; i<nMC;i++){
	   if(abs((*mcPID)[i]) == PID){
	     if(!( mcStatus->at(i) == 62)) continue;
	     Wfound=true;
	     bosonPID = (*mcPID)[i];
	     bosonPt = (*mcPt)[i];
	   }
	 }
	 if ( bosonPt > 0 ){
	   if(debug)cout<<"Accessing nlo ewk, qcd"<<endl;
	   nlo_ewk = NLO_EWK->GetBinContent(NLO_EWK->GetXaxis()->FindBin(bosonPt));
	   nlo_qcd_binned=NLO_QCD->GetBinContent(NLO_QCD->GetXaxis()->FindBin(bosonPt));
	   if (found_Wjet_sample) {
	     nlo_qcd = exponential(bosonPt,1.053, 3.163e-3, 0.746);
	   } else if (found_DYjet_sample) {
	     //nlo_qcd = exponential(bosonPt,1.434, 2.210e-3, 0.443);
	     nlo_qcd = NLO_QCD->GetBinContent(NLO_QCD->GetXaxis()->FindBin(bosonPt));
	   } //else if (type == GJets) {
	   //nlo_qcd = exponential(bosonPt,1.159, 1.944e-3, 1.0);
	   // }
	   
	   if(debug)cout<<"Accessing nnlo qcd"<<endl;
	   nnlo_qcd = NNLO_QCD->GetBinContent(NNLO_QCD->GetXaxis()->FindBin(bosonPt));
	 }
	 // if (isNLO) kfactor = nlo_ewk * nnlo_qcd;
	 // else kfactor = nlo_ewk * nlo_qcd * nnlo_qcd;
      	 if(debug) cout<<"apply kfactor"<<endl;
	 if(nlo_ewk*nlo_qcd !=0 ) kfactor = nlo_ewk * nlo_qcd;
	 
       }
       weight=weight*kfactor;
       //cout<<" kfactor = "<<kfactor<<endl;
       event_weight=weight;
       if(debug)cout<<"reco selections begin"<<endl;
       
       muCand.clear(); mu2Cand.clear();  tauCand.clear();
       event_weight=weight;
       ////// DY Z-> ll signal region -  isolated begin
       if(muMu_selector)
	 {
	   if(metFilters==0 )
	     {
	       if(debug)cout<<"metfilters selected"<<endl;
	       if(isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	       nMETFiltersPassed_dyll+=event_weight;
	       if(debug)cout<<"genweight applied"<<endl;
	       if(   passSingleTriggerPaths || passCrossTrigger )
		 {
		   nSingleTrgPassed_dyll+=event_weight;
		   if(debug)cout<<"trigger selected"<<endl;
		   muCand = getMuCand(20,2.1);  ///// muons selected 
		   if( muCand.size() >0 ) 
		     { 
		       nGoodMuonPassed_dyll+=event_weight;
		       if(debug)cout<<"this worked Line 443"<<endl;
		       std::vector<int> iMuPlus; iMuPlus.clear(); 
		       std::vector<int> iMuMinus; iMuMinus.clear();
		       for(int i=0; i<muCand.size(); i++){
			 if(muCharge->at(muCand[i]) < 0) iMuMinus.push_back(muCand[i]);
			 if(muCharge->at(muCand[i]) > 0) iMuPlus.push_back(muCand[i]);
		       }
		       if(iMuPlus.size()>0 && iMuMinus.size()>0)
			 {
			   //mu2Cand = getMu2Cand(20,2.4, muCand[0]);
			   //if( mu2Cand.size() >0 )
			   //{
			   //if( muCharge->at(muCand[0])*muCharge->at(mu2Cand[0]) <0 )
			   {
			     tauCand = getTauCand(30,2.3);
			     if( tauCand.size()>0 ) 
			       {
				 nGoodTauPassed_dyll+=event_weight;
				 if(debug)cout<<"this worked Line 305"<<endl;
				 reco_mu.clear();reco_tau.clear(); reco_mu2.clear();
				 if( muPt->at(iMuMinus[0]) > muPt->at(iMuPlus[0])   ) { reco_mu=iMuMinus; reco_mu2=iMuPlus; }
				 else { reco_mu=iMuPlus; reco_mu2=iMuMinus;}
				 reco_tau=tauCand;
				 //reco_mu=muCand; reco_tau=tauCand; reco_mu2=mu2Cand;
				 if ( MatchTriggerFilter(reco_mu[0], reco_tau[0]) )
				   {
				     if ( muCharge->at(reco_mu[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
				       {
					 nGoodMuTauPassed_dyll+=event_weight;
					 muonP4.SetPtEtaPhiE(muPt->at(reco_mu[0]), muEta->at(reco_mu[0]), muPhi->at(reco_mu[0]), muE->at(reco_mu[0]));
					 tauP4.SetPtEtaPhiE(tau_Pt->at(reco_tau[0]), tau_Eta->at(reco_tau[0]), tau_Phi->at(reco_tau[0]), tau_Energy->at(reco_tau[0]));
					 if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
					 
					 muonP4 = muonP4*muRC_sf;
					 if(debug)cout<<" sf : "<<getScaleFactors( reco_mu[0] , reco_tau[0] , false , isMC , debug ) <<endl;
					 if (isMC) event_weight = event_weight * getScaleFactors( reco_mu[0] , reco_tau[0] , false , isMC , debug );
					 afterSF4+=event_weight;
					 if( thirdLeptonVeto() < 0 )
					   {
					     nPassedThirdLepVeto_dyll+=event_weight;
					     
					     if( passBjetVeto() == true)
					       {
						 nPassedBjetVeto_dyll+=event_weight;
						 
						 double deltaR = delta_R( muPhi->at(reco_mu[0]),  muEta->at(reco_mu[0]),  muPhi->at(reco_mu2[0]),  muEta->at(reco_mu2[0]));
						 if(deltaR > 0.5 )
						   {
						     nDeltaRPassed_dyll+=event_weight;
						     if(isMC==false)event_weight=1.0;
						     if(debug)cout<<"this worked Line 374"<<endl;
						     fillHist_dyll("5_dyll",  reco_mu[0], reco_mu2[0], reco_tau[0], event_weight);
						     double mT_muMet = TMass_F((muPt->at(reco_mu[0])),(muPhi->at(reco_mu[0])),pfMET,pfMETPhi  );
						     if( mT_muMet < 50)
						       {
							 fillHist_dyll("6_dyll", reco_mu[0], reco_mu2[0], reco_tau[0], event_weight);
						       }
						   }
					       }
					   }
				       }
				   }
			       }
			   }
			 }
		     }
		 }
	     }  //////// signal region DY Z->ll end
       ////// fake background region - antiisolated  DY Z->ll begin
       if(debug)cout<<"moving to fake bkg"<<endl;
       event_weight=weight;
       muCand.clear(); tauCand.clear();  mu2Cand.clear();
       if(metFilters==0 )
	 {
	   if(isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed_dyll_fr+=event_weight;
	   if(  passSingleTriggerPaths || passCrossTrigger  )
	     {
	       nSingleTrgPassed_dyll_fr+=event_weight;
	       if(debug)cout<<"trigger selected"<<endl;
	       muCand = getMuCand(20,2.1);  ///// muons selected 
	       if( muCand.size() >0 ) 
		 { 
		   nGoodMuonPassed_dyll_fr+=event_weight;
		   if(debug)cout<<"this worked Line 514"<<endl;
		   std::vector<int> iMuPlus; iMuPlus.clear();
		   std::vector<int> iMuMinus; iMuMinus.clear();
                   for(int i=0; i<muCand.size(); i++){
                     if(muCharge->at(muCand[i]) < 0) iMuMinus.push_back(muCand[i]);
                     if(muCharge->at(muCand[i]) > 0) iMuPlus.push_back(muCand[i]);
                   }
                   if(iMuPlus.size()>0 && iMuMinus.size()>0)
                     {
		       //mu2Cand = getMu2Cand(20,2.4, muCand[0]);
		       //if( mu2Cand.size() >0 )
		       //{
                       //if( muCharge->at(muCand[0])*muCharge->at(mu2Cand[0]) <0  )
		       {
			 tauCand = getAISRTauCand(30,2.3);
			 if( tauCand.size()>0 ) 
			   {
			     nGoodTauPassed_dyll_fr+=event_weight;
			     if(debug)cout<<"fr tau selection passed"<<endl;
			     reco_mu.clear();reco_tau.clear(); reco_mu2.clear();
			     if( muPt->at(iMuMinus[0]) > muPt->at(iMuPlus[0])   ) { reco_mu=iMuMinus; reco_mu2=iMuPlus; }
			     else { reco_mu=iMuPlus; reco_mu2=iMuMinus;}
			     reco_tau=tauCand;
			     //reco_mu=muCand; reco_tau=tauCand; reco_mu2=mu2Cand;
			     if ( MatchTriggerFilter(reco_mu[0], reco_tau[0]) )
			       {
				 if ( muCharge->at(reco_mu[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
				   {
				     nGoodMuTauPassed_dyll_fr+=event_weight;
				     muonP4.SetPtEtaPhiE(muPt->at(reco_mu[0]), muEta->at(reco_mu[0]), muPhi->at(reco_mu[0]), muE->at(reco_mu[0]));
				     tauP4.SetPtEtaPhiE(tau_Pt->at(reco_tau[0]), tau_Eta->at(reco_tau[0]), tau_Phi->at(reco_tau[0]), tau_Energy->at(reco_tau[0]));
				     
				     muonP4 = muonP4*muRC_sf;
				     event_weight = event_weight* getFR(reco_tau[0]);
				     
				     /////
				     if(debug)cout<<" fake bkg sf : "<<getScaleFactors( reco_mu[0] , reco_tau[0] , true , isMC , debug ) <<endl;
				     if(isMC) event_weight = event_weight * getScaleFactors( reco_mu[0] , reco_tau[0] , true , isMC , debug );
				     /////
				     if( thirdLeptonVeto() < 0 )
				       {
					 nPassedThirdLepVeto_dyll_fr+=event_weight;
					 
					 if( passBjetVeto() == true)
					   {
					     nPassedBjetVeto_dyll_fr+=event_weight;
					     
					     double deltaR = delta_R( muPhi->at(reco_mu[0]),  muEta->at(reco_mu[0]),  muPhi->at(reco_mu2[0]),  muEta->at(reco_mu2[0]));
					     if(deltaR > 0.5 )
					       {
						 nDeltaRPassed_dyll_fr+=event_weight;
						 if(debug)cout<<"this worked Line 425"<<endl;
						 if(debug)cout<<"evnt_weight = "<<event_weight<<endl;
						 fillHist_dyll("5_dyll_fr", reco_mu[0], reco_mu2[0], reco_tau[0], event_weight);
						 double mT_muMet = TMass_F((muPt->at(reco_mu[0])),(muPhi->at(reco_mu[0])),pfMET,pfMETPhi  );
						 if( mT_muMet < 50 ) 
						   {
						     fillHist_dyll("6_dyll_fr", reco_mu[0], reco_mu2[0], reco_tau[0], event_weight);
						   }
					       }
					   }
				       }
				   }
			       }
			   }
		       }
		     }
		 }
	     }
	 }
	 }
       ////// fake rate anti isolated region Dy Z->ll  end

       ////// signal region -  isolated begin
       event_weight=weight;
       muCand.clear(); tauCand.clear();  mu2Cand.clear();
       if(muTau_selector)
	 {
       if(metFilters==0)
	 {
	   
	   if(debug)cout<<"metfilters selected"<<endl;
	   if(isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed+=event_weight;
	   if(debug)cout<<"genweight applied"<<endl;
	   if(  passSingleTriggerPaths || passCrossTrigger )
	     {
	       nSingleTrgPassed+=event_weight;
	       if(debug)cout<<"trigger selected"<<endl;
	       muCand = getMuCand(20,2.1);  ///// muons selected 
	       if( muCand.size() >0 ) 
		 { 
		   nGoodMuonPassed+=event_weight;
		   if(debug)cout<<"this worked Line 284"<<endl;
		   tauCand = getTauCand(30,2.3);
		   if( tauCand.size()>0 ) 
		     {
		       nGoodTauPassed+=event_weight;
		       if(debug)cout<<"this worked Line 305"<<endl;
		       
		       reco_mu.clear();reco_tau.clear();
		       reco_mu=muCand; reco_tau=tauCand;
		       
		       if ( MatchTriggerFilter(reco_mu[0], reco_tau[0]) )
			 {
			   if ( muCharge->at(reco_mu[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
			     {
			       nGoodMuTauPassed+=event_weight;
			       muonP4.SetPtEtaPhiE(muPt->at(reco_mu[0]), muEta->at(reco_mu[0]), muPhi->at(reco_mu[0]), muE->at(reco_mu[0]));
			       tauP4.SetPtEtaPhiE(tau_Pt->at(reco_tau[0]), tau_Eta->at(reco_tau[0]), tau_Phi->at(reco_tau[0]), tau_Energy->at(reco_tau[0]));
			       if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
			       
			       muonP4 = muonP4*muRC_sf;
			       if(debug)cout<<" sf : "<<getScaleFactors( reco_mu[0] , reco_tau[0] , false , isMC , debug ) <<endl;
			       if (isMC) event_weight = event_weight * getScaleFactors( reco_mu[0] , reco_tau[0] , false , isMC , debug );
			       afterSF4+=event_weight;
			       if( thirdLeptonVeto() < 0 )
				 {
				   nPassedThirdLepVeto+=event_weight;
				   
				   if( passBjetVeto() == true)
				     {
				       nPassedBjetVeto+=event_weight;
				       
				       double deltaR = delta_R(muonP4.Phi(), muonP4.Eta(), tauP4.Phi(), tauP4.Eta());
				       if(deltaR > 0.5 )
					 {
					   nDeltaRPassed+=event_weight;
					   if(isMC==false)event_weight=1.0;
					   if(debug)cout<<"this worked Line 374"<<endl;
					   fillHist("5", muonP4, tauP4, reco_mu[0], reco_tau[0], event_weight);
					   double mT_muMet = TMass_F((muPt->at(reco_mu[0])),(muPhi->at(reco_mu[0])),pfMET,pfMETPhi  );
					   if( mT_muMet < 50)
					     {
					       fillHist("6", muonP4, tauP4, reco_mu[0], reco_tau[0], event_weight);
					     }
					 }
				     }
				 }
			     }
			 }
		     }
		 }
	     }
	 }  //////// signal region end
	 
       ////// fake background region - antiisolated begin
       if(debug)cout<<"moving to fake bkg"<<endl;
       event_weight=weight;
       muCand.clear(); tauCand.clear();
       if(metFilters==0)
	 {
	   
	   if(isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed_fr+=event_weight;
	   if(  passSingleTriggerPaths || passCrossTrigger  )
	     {
	       nSingleTrgPassed_fr+=event_weight;
	       if(debug)cout<<"trigger selected"<<endl;
	       muCand = getMuCand(20,2.1);  ///// muons selected 
	       if( muCand.size() >0 ) 
		 { 
		   nGoodMuonPassed_fr+=event_weight;
		   if(debug)cout<<"this worked Line 373"<<endl;
		   tauCand = getAISRTauCand(30,2.3);
		   if( tauCand.size()>0 ) 
		     {
		       nGoodTauPassed_fr+=event_weight;
		       if(debug)cout<<"fr tau selection passed"<<endl;
		       
		       reco_mu.clear();reco_tau.clear();
		       reco_mu=muCand; reco_tau=tauCand;
		       if ( MatchTriggerFilter(reco_mu[0], reco_tau[0]) )
			 {
			   
			   if ( muCharge->at(reco_mu[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
			     {
			       nGoodMuTauPassed_fr+=event_weight;
			       muonP4.SetPtEtaPhiE(muPt->at(reco_mu[0]), muEta->at(reco_mu[0]), muPhi->at(reco_mu[0]), muE->at(reco_mu[0]));
			       tauP4.SetPtEtaPhiE(tau_Pt->at(reco_tau[0]), tau_Eta->at(reco_tau[0]), tau_Phi->at(reco_tau[0]), tau_Energy->at(reco_tau[0]));
			       
			       
			       muonP4 = muonP4*muRC_sf;
			       
			       event_weight = event_weight* getFR(reco_tau[0]);
			       
			       /////
			       if(debug)cout<<" fake bkg sf : "<<getScaleFactors( reco_mu[0] , reco_tau[0] , true , isMC , debug ) <<endl;
			       if(isMC) event_weight = event_weight * getScaleFactors( reco_mu[0] , reco_tau[0] , true , isMC , debug );
			       /////
			       if( thirdLeptonVeto() < 0 )
				 {
				   nPassedThirdLepVeto_fr+=event_weight;
				   
				   if( passBjetVeto() == true)
				     {
				       nPassedBjetVeto_fr+=event_weight;
				       
				       double deltaR = delta_R(muonP4.Phi(), muonP4.Eta(), tauP4.Phi(), tauP4.Eta());
				       if(deltaR > 0.5 )
					 {
					   nDeltaRPassed_fr+=event_weight;
					   if(debug)cout<<"this worked Line 425"<<endl;
					   if(debug)cout<<"evnt_weight = "<<event_weight<<endl;
					   fillHist("5_fr", muonP4, tauP4, reco_mu[0], reco_tau[0], event_weight);
					   double mT_muMet = TMass_F((muPt->at(reco_mu[0])),(muPhi->at(reco_mu[0])),pfMET,pfMETPhi  );
					   if( mT_muMet < 50 ) 
					     {
					       fillHist("6_fr", muonP4, tauP4, reco_mu[0], reco_tau[0], event_weight);
					     }
					 }
				     }
				 }
			     }
			 }
		     }
		 }
	     }
	 }
	 }
       ////// fake rate anti isolated region end

       
       report_test = nentriesToCheck/20;
       while (report_test>10)
	 {
	   report_test=report_test/10;
	   report_i++;
	 }
       if(nentriesToCheck>20)
	 reportEvery = report_test*pow(10,report_i);
       else 
	 reportEvery = 1;
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
   

     
   h_cutflow_n->SetBinContent(1,nInspected_genWeighted );
   h_cutflow_n->SetBinContent(2, nSingleTrgPassed);
   h_cutflow_n->SetBinContent(3, nGoodMuonPassed);
   h_cutflow_n->SetBinContent(4, nGoodTauPassed);
   h_cutflow_n->SetBinContent(5, nGoodMuTauPassed);
   h_cutflow_n->SetBinContent(6, nPassedThirdLepVeto);
   h_cutflow_n->SetBinContent(7, nPassedBjetVeto);
   h_cutflow_n->SetBinContent(8, nDeltaRPassed);
   
   h_cutflow_n_fr->SetBinContent(1,nInspected_genWeighted );
   h_cutflow_n_fr->SetBinContent(2, nSingleTrgPassed_fr);
   h_cutflow_n_fr->SetBinContent(3, nGoodMuonPassed_fr);
   h_cutflow_n_fr->SetBinContent(4, nGoodTauPassed_fr);
   h_cutflow_n_fr->SetBinContent(5, nGoodMuTauPassed_fr);
   h_cutflow_n_fr->SetBinContent(6, nPassedThirdLepVeto_fr);
   h_cutflow_n_fr->SetBinContent(7, nPassedBjetVeto_fr);
   h_cutflow_n_fr->SetBinContent(8, nDeltaRPassed_fr);

   /// dy Z->ll
   h_cutflow_n_dyll->SetBinContent(1,nInspected_genWeighted );
   h_cutflow_n_dyll->SetBinContent(2, nSingleTrgPassed_dyll);
   h_cutflow_n_dyll->SetBinContent(3, nGoodMuonPassed_dyll);
   h_cutflow_n_dyll->SetBinContent(4, nGoodTauPassed_dyll);
   h_cutflow_n_dyll->SetBinContent(5, nGoodMuTauPassed_dyll);
   h_cutflow_n_dyll->SetBinContent(6, nPassedThirdLepVeto_dyll);
   h_cutflow_n_dyll->SetBinContent(7, nPassedBjetVeto_dyll);
   h_cutflow_n_dyll->SetBinContent(8, nDeltaRPassed_dyll);
   
   h_cutflow_n_dyll_fr->SetBinContent(1,nInspected_genWeighted );
   h_cutflow_n_dyll_fr->SetBinContent(2, nSingleTrgPassed_dyll_fr);
   h_cutflow_n_dyll_fr->SetBinContent(3, nGoodMuonPassed_dyll_fr);
   h_cutflow_n_dyll_fr->SetBinContent(4, nGoodTauPassed_dyll_fr);
   h_cutflow_n_dyll_fr->SetBinContent(5, nGoodMuTauPassed_dyll_fr);
   h_cutflow_n_dyll_fr->SetBinContent(6, nPassedThirdLepVeto_dyll_fr);
   h_cutflow_n_dyll_fr->SetBinContent(7, nPassedBjetVeto_dyll_fr);
   h_cutflow_n_dyll_fr->SetBinContent(8, nDeltaRPassed_dyll_fr);
   ///
   
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
  //cout<<"cloning nEvents hist"<<endl;
  h_nEvents = (TH1F*)((TH1F*)file_in->Get("nEvents"))->Clone(TString("nEvents"));
  file_in->Close();
  

}

//Fill the sequential histos at a particular spot in the sequence


void mutau_analyzer::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{
  
  //  h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);

}


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
       bool filter = false;
       if( muPt->at(iMu) > muPtCut  && fabs(muEta->at(iMu))< muEtaCut  && fabs(muDz->at(iMu)) < 0.2 && fabs(muD0->at(iMu))<0.045 ) kinematic = true;
       if( muIDbit->at(iMu)>>1&1==1 ) muonID = true;//|| muIDbit->at(iMu)>>8&1==1 || muIDbit->at(iMu)>>17&1==1  ) muonID = true;
       float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
       if( relMuIso < 0.15 ) muonIso = true;
       
       if( ( HLTEleMuX>>21&1 == 1 && muPt->at(iMu) > 28 ) 
	   || ( HLTEleMuX>>60&1 == 1 && muPt->at(iMu) > 25 )
	   || ( HLTTau>>0&1 == 1 && muPt->at(iMu) <= 25  && muPt->at(iMu) > 21 && fabs(muEta->at(iMu))< 2.1)
	   //|| ( HLTTau>>8&1 == 1 && muPt->at(iMu) <= 25  && muPt->at(iMu) > 21 && fabs(muEta->at(iMu))< 2.1)
	   ) trigger = true;
       if( kinematic==true  &&  muonID==true &&  muonIso==true && trigger==true ){
	 tmpCand.push_back(iMu);
       }                                                                                      
     }                                                                                       
  return tmpCand;
  
}

std::vector<int> mutau_analyzer::getMu2Cand(double muPtCut, double muEtaCut, int mu1Index){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over muons                                                                     
   for(int iMu=0;iMu<nMu;iMu++)
     {
       if(iMu==mu1Index)
	 continue;
       bool kinematic = false;
       bool muonID = false;
       bool muonIso =  false;
       bool trigger = false;
       bool filter = false;
       if( muPt->at(iMu) > muPtCut  && fabs(muEta->at(iMu))< muEtaCut  && fabs(muDz->at(iMu)) < 0.2 && fabs(muD0->at(iMu))<0.045 ) kinematic = true;
       if( muIDbit->at(iMu)>>1&1==1 ) muonID = true;//|| muIDbit->at(iMu)>>8&1==1 || muIDbit->at(iMu)>>17&1==1  ) muonID = true;
       float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
       if( relMuIso < 0.15 ) muonIso = true;
       if( ( HLTEleMuX>>21&1 == 1 && muPt->at(iMu) > 28 )
	   || ( HLTEleMuX>>60&1 == 1 && muPt->at(iMu) > 25  )
	   || ( HLTTau>>0&1 == 1 && muPt->at(iMu) <= 25  && muPt->at(iMu) > 21 && fabs(muEta->at(iMu))< 2.1)
	   //|| ( HLTTau>>8&1 == 1 && muPt->at(iMu) <= 25  && muPt->at(iMu) > 21 && fabs(muEta->at(iMu))< 2.1)
	   ) trigger = true;
       
       if( kinematic==true  &&  muonID==true &&  muonIso==true && trigger==true ){
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
      bool trigger = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool mutau_separation=false;
      bool newDecayModeFinding=false;
      bool tau_reject=false;
      bool filter = false;
      if( tau_Pt->at(iTau) > tauPtCut 
	  && fabs( tau_Eta->at(iTau))< tauEtaCut
	  && tau_LeadChargedHadron_dz->at(iTau) < 0.2
  	  )kinematic = true;
      if( tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1 ) tauIsolation=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byTightDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;

      if( ( HLTEleMuX>>21&1 == 1 )
	  || ( HLTEleMuX>>60&1 == 1 )
	  || ( HLTTau>>0&1 == 1 && tau_Pt->at(iTau) >32 && abs(tau_Eta->at(iTau)) < 2.1 )
	  //|| ( HLTTau>>8&1 == 1 && tau_Pt->at(iTau) >32 && abs(tau_Eta->at(iTau)) < 2.1 )
	  ) trigger=true;
      if( kinematic==true    
	  && decayModeCut==true   
	  && tauIsolation==true 
	  && tau_reject==true   
	  && newDecayModeFinding==true
	  && trigger==true
	  )
	{
	  tmpCand.push_back(iTau);
    	}                                                           
    }                                                                                       
  return tmpCand;
  
}
std::vector<int> mutau_analyzer::getAISRTauCand(double tauPtCut, double tauEtaCut){
  std::vector<int> tmpCand;  tmpCand.clear();
  for(int iTau=0;iTau<nTau;iTau++) //Loop over taus
    {
      bool kinematic = false;
      bool trigger = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool mutau_separation=false;
      bool newDecayModeFinding=false;
      bool tau_reject=false;
      bool filter = false;
      if( tau_Pt->at(iTau) > tauPtCut 
	  && fabs( tau_Eta->at(iTau))< tauEtaCut 
	  && tau_LeadChargedHadron_dz->at(iTau) < 0.2
  	  )kinematic = true;
      if( tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 && tau_byMediumDeepTau2017v2p1VSjet->at(iTau)!=1 ) tauIsolation=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byTightDeepTau2017v2p1VSmu->at(iTau)==1 )tau_reject=true;
      
      if( ( HLTEleMuX>>21&1 == 1 )
	  || ( HLTEleMuX>>60&1 == 1 )
	  || ( HLTTau>>0&1 == 1 && tau_Pt->at(iTau) >32 && abs(tau_Eta->at(iTau)) < 2.1 )
	  //|| ( HLTTau>>8&1 == 1 && tau_Pt->at(iTau) >32 && abs(tau_Eta->at(iTau)) < 2.1 )
	  ) trigger=true;

      if( kinematic==true    
	  && decayModeCut==true   
	  && tauIsolation==true 
	  && tau_reject==true   
	  && newDecayModeFinding==true
	   && trigger==true
     	  )
	{
	  tmpCand.push_back(iTau);
    	}                                                           
    }                                                                                       
  return tmpCand;  
}
std::vector<int> mutau_analyzer::getJetCand(int muIndex, int tauIndex){
  std::vector<int> tmpCand;  tmpCand.clear();
  for(int iJet=0;iJet<nJet;iJet++) //Loop over jets
    {
      bool kinematic30 = false;
      bool kinematic50 = false; bool kinematic50Loose = false;
      bool jet_id = false; bool drPassed=false;
      if( jetPt->at(iJet) > 50 
	  && abs(jetEta->at(iJet))<2.65
	  && abs(jetEta->at(iJet))>3.139
	  && (jetID->at(iJet)>>0&1)==1
	  ) kinematic50=true;
      // else if( jetPt->at(iJet) < 50  
      // 	       && jetPUFullID->at(iJet)>>1&1==1 
      // 	       ) kinematic50Loose=true;
      else if( jetPt->at(iJet) > 30
	       && abs(jetEta->at(iJet))<4.7
	       && (jetID->at(iJet)>>0&1)==1
	       ) kinematic30=true;
                  
      double dr_jetMu=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , muPhi->at(muIndex), muEta->at(muIndex) );
      double dr_jetTau=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , tau_Phi->at(tauIndex), tau_Eta->at(tauIndex) );
      if( dr_jetMu>0.5 && dr_jetTau>0.5 )
	drPassed=true;
	  
      if( (kinematic50 || kinematic30 ) && drPassed==true)
	tmpCand.push_back(iJet);
    }
  return tmpCand;
}
//The noisy jets are defined as: 20 < pt < 50 && abs(eta) > 2.65 && abs(eta) < 3.139. 
bool mutau_analyzer::noisyJet2017(){

  bool noisyJet = false;
  for(int iJet=0;iJet<nJet;iJet++) //Loop over jets
    {
      bool kinematic = false;
      if( jetRawPt->at(iJet) > 20 
	  && jetRawPt->at(iJet) < 50
	  && abs(jetEta->at(iJet)) > 2.65
          && abs(jetEta->at(iJet)) < 3.139
	  ) kinematic=true;

      if( kinematic )
        noisyJet=true;
    }
  return noisyJet;
}
int mutau_analyzer::thirdLeptonVeto(){
  std::vector<int> tmpCand;
  tmpCand.clear();
  int thirdLepIndex = -1;
  bool thirdLepVeto=true;
  for(int iEle=0; iEle < nEle;iEle++)
    {
      bool kinematic = false;
      if( (*elePt)[iEle] > 10.0  
	  && fabs((*eleEta)[iEle])< 2.5 
	  && (*eleD0)[iEle] < 0.045 
	  && (*eleDz)[iEle] < 0.2 
	  && eleMissHits->at(iEle) <= 1 && eleConvVeto->at(iEle)==1
	  ) kinematic = true;
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
  return  sqrt(2.0*LepPt*met*(1.0-cos(DeltaPhi(LepPhi, metPhi))));
  //return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
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
float mutau_analyzer::pTvecsum_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float pt_vecSum = (a + b+ met).Pt();
  return pt_vecSum;
}

bool mutau_analyzer::passBjetVeto()
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  // for(int iJets=0; iJets<nJet ; iJets++){
  //   if(jetCSV2BJetTags->at(iJets) > 0.8838) tmpCand++;
  //   // CSV B jet tag for selecting bJets is medium WP (jetCSV2BJetTags > 0.8838.)
  // }
  for(int iJets=0; iJets<nJet ; iJets++){
    if( jetPt->at(iJets) > 25  && abs(jetEta->at(iJets)) < 2.4 && jetID->at(iJets)>>0&1==1 ){
      tmpJetCand.push_back(iJets);
    }
  }
  if(tmpJetCand.size() >=1 ){
    // atleast one jet ==> events pass medium 
    if( (jetDeepCSVTags_b->at(tmpJetCand[0]) + jetDeepCSVTags_bb->at(tmpJetCand[0])) > 0.4941  )
      return veto = false;
  }
  else if(tmpJetCand.size() >= 2){
    // atleast 2 jets ==> events pass loose
    if( (jetDeepCSVTags_b->at(tmpJetCand[0]) + jetDeepCSVTags_bb->at(tmpJetCand[0])) > 0.1522   )
      return veto = false;
  }
  return veto;
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
std::vector<int> mutau_analyzer::getGenMu(){
  std::vector<int> tmpCand;
  tmpCand.clear();
  int count1=0; int count2=0;
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>2&1==1 || genMatch2->at(imc)>>4&1==1 ) { count1++;}
    if( (genMatch2->at(imc)>>2&1==1 || genMatch2->at(imc)>>4&1==1) && fabs(mcPID->at(imc))==13 ){  tmpCand.push_back(imc); count2++;}
  }
  //cout<<"count1:"<<count1<<"  count2:"<<count2<<endl;
  return tmpCand; 
}
bool mutau_analyzer::hasGenTau(){
  bool found_genTau=false;
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>5&1==1) {  found_genTau=true;}
  }
  return found_genTau;
}
float mutau_analyzer::exponential(float x,float a,float b,float c) {
  return a * TMath::Exp(-b * x) + c;
}
double mutau_analyzer::getScaleFactors( int muIndex , int tauIndex, bool fakeBkg , bool isMC, bool debug)
{
  double rv_sf=1.0;
  double sf_IsoEff = 1.0; 
  double sf_muTrg = 1.0;
  double sf_muID = 1.0;
  double sf_tauidSF_m = 1.0;
  double sf_tauidSF_vvvl = 1.0;
  double sf_tauesSF = 1.0;
  double sf_fakeEle = 1.0; double sf_fakeMu = 1.0;
  double sf_taufesSF = 1.0;
  int genMatchTau = 0;
  if(isMC) genMatchTau = gen_matching();
  double recoMuonPt=0.0;
  if (muPt->at(muIndex) < 120)
    recoMuonPt=muPt->at(muIndex);
  else
    recoMuonPt = 119;
  sf_muID = h_muIDSF->GetBinContent(h_muIDSF->GetXaxis()->FindBin(recoMuonPt),h_muIDSF->GetYaxis()->FindBin(abs(muEta->at(muIndex))));
  sf_IsoEff = h_muIsoSF->GetBinContent(h_muIDSF->GetXaxis()->FindBin(recoMuonPt),h_muIDSF->GetYaxis()->FindBin(abs(muEta->at(muIndex))));  
  sf_muTrg = h_muTrgSF->GetBinContent(h_muIDSF->GetXaxis()->FindBin(muPt->at(muIndex)),h_muIDSF->GetYaxis()->FindBin(abs(muEta->at(muIndex))));

  sf_tauidSF_m = h_tauidSF_m->GetBinContent(h_tauidSF_m->GetXaxis()->FindBin(tau_DecayMode->at(tauIndex)));
  //sf_tauidSF_m = fn_tauIDSF_m->Eval(tau_Pt->at(reco_tau[0]));
  sf_tauidSF_vvvl = h_tauidSF_vvvl->GetBinContent(h_tauidSF_vvvl->GetXaxis()->FindBin(tau_DecayMode->at(tauIndex)));
  //sf_tauidSF_vvvl = fn_tauIDSF_vvl->Eval(tau_Pt->at(reco_tau[0]));
  sf_tauesSF = h_tauesSF->GetBinContent(h_tauesSF->GetXaxis()->FindBin(tau_DecayMode->at(tauIndex)));

  //if(genMatchTau==1 || genMatchTau==3)
  //sf_fakeEle = h_tauFakeEleSF->GetBinContent(h_tauFakeEleSF->GetXaxis()->FindBin(abs(tau_Eta->at(tauIndex))));
  //if(genMatchTau==2 || genMatchTau==4)
  //sf_fakeMu = h_tauFakeMuSF->GetBinContent(h_tauFakeMuSF->GetXaxis()->FindBin(abs(tau_Eta->at(tauIndex))));
  if(genMatchTau==2 || genMatchTau==4){
    if(tau_DecayMode->at(tauIndex)==0)
      {
	if(abs(tau_Eta->at(tauIndex)) < 0.4 ) sf_fakeMu=1.14;
	if(abs(tau_Eta->at(tauIndex)) > 0.4 
	   && abs(tau_Eta->at(tauIndex)) < 0.8 ) sf_fakeMu=1.0;
	if(abs(tau_Eta->at(tauIndex)) > 0.8
           && abs(tau_Eta->at(tauIndex)) < 1.2 ) sf_fakeMu=0.87;
	if(abs(tau_Eta->at(tauIndex)) > 1.2
           && abs(tau_Eta->at(tauIndex)) < 1.7 ) sf_fakeMu=0.52;
	if(abs(tau_Eta->at(tauIndex)) > 1.7
           && abs(tau_Eta->at(tauIndex)) < 2.3 ) sf_fakeMu=1.47;
      }
    if(tau_DecayMode->at(tauIndex)==1)
      {
	if(abs(tau_Eta->at(tauIndex)) > 0.0
           && abs(tau_Eta->at(tauIndex)) < 0.4 ) sf_fakeMu=0.69;
      }
  }
  if(genMatchTau==1 || genMatchTau==3){
    if(tau_DecayMode->at(tauIndex)==0)
      {
	if(abs(tau_Eta->at(tauIndex)) < 1.479 ) sf_fakeEle=0.98;
	if(abs(tau_Eta->at(tauIndex)) > 1.479 ) sf_fakeEle=0.80;
      }
    if(tau_DecayMode->at(tauIndex)==1)
      {
	if(abs(tau_Eta->at(tauIndex)) < 1.479 ) sf_fakeEle=1.07;
        if(abs(tau_Eta->at(tauIndex)) > 1.479 ) sf_fakeEle=0.64;
      }
  }
  if(tau_DecayMode->at(tauIndex)==0 && abs(tau_Eta->at(tauIndex))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(1);
  if(tau_DecayMode->at(tauIndex)==0 && abs(tau_Eta->at(tauIndex))>1.4 )  sf_taufesSF = h_taufesSF->Eval(3);
  if(tau_DecayMode->at(tauIndex)==1 && abs(tau_Eta->at(tauIndex))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(5);
  if(tau_DecayMode->at(tauIndex)==1 && abs(tau_Eta->at(tauIndex))>1.4 )  sf_taufesSF = h_taufesSF->Eval(7);
  
  
  //event_weight=event_weight * sf_muID * sf_IsoEff * sf_muTrg * sf_tauidSF_m * sf_tauesSF * sf_fakeEle * (sf_fakeMu);
  // if(fakeBkg)
  //   rv_sf = sf_muID * sf_IsoEff * sf_muTrg * sf_tauidSF_m * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF ;
  // else
  //   rv_sf = sf_muID * sf_IsoEff * sf_muTrg * sf_tauidSF_m * sf_tauidSF_vvvl * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF ;
  if(fakeBkg)
    rv_sf = sf_muID * sf_IsoEff * sf_muTrg * sf_tauidSF_m * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF;
  else
    rv_sf = sf_muID * sf_IsoEff * sf_muTrg * sf_tauidSF_m * sf_tauidSF_vvvl * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF;
  
  return rv_sf;

}

bool mutau_analyzer::MatchTriggerFilter(int muIndex, int tauIndex)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool passFilter = true;
  bool muTriggerFilterMatch=false;
  int nMuTriggerFilterMatch=0;
  bool tauTriggerFilterMatch=false;
  int nTauTriggerFilterMatch=0;
  for(int ifilter=33;ifilter<56;ifilter++)
    {
      if(muFiredTrgs->at(muIndex)>>ifilter&1==1)
	{
	  muTriggerFilterMatch=true;
	  nMuTriggerFilterMatch++;
	}
    }
  for(int ifilter=0;ifilter<18;ifilter++)
    {
      if(tauFiredTrgs->at(tauIndex)>>ifilter&1==1)
	{
	  tauTriggerFilterMatch=true;
	  nTauTriggerFilterMatch++;
	}
    }
  if( (HLTEleMuX>>21&1 == 1 && nMuTriggerFilterMatch==23 )
      || (HLTEleMuX>>60&1 == 1 && nMuTriggerFilterMatch==23 ) 
      || (HLTTau>>0&1 ==1  && nMuTriggerFilterMatch==23 )
      )
    passFilter=true;
    
  return passFilter;
}



double  mutau_analyzer::getFR(int tauIndex){
  double frWeight=1.0;
  double tau_FR = 1.0;
  double tauPt=0.0;
  if( tau_Pt->at(tauIndex) < 120 )
    tauPt=tau_Pt->at(tauIndex);
  else
    tauPt=119.0;
  if ( tau_DecayMode->at(tauIndex)==0 )
    {
      tau_FR = h_tauFR_0->Eval(tauPt);
      frWeight = tau_FR/(1-tau_FR);
    }
  
  if ( tau_DecayMode->at(tauIndex)==1 )
    {
      tau_FR = h_tauFR_1->Eval(tauPt);
      frWeight = tau_FR/(1-tau_FR);
    }
  
  if ( tau_DecayMode->at(tauIndex)==10 )
    {
      tau_FR = h_tauFR_10->Eval(tauPt);
      frWeight = tau_FR/(1-tau_FR);
    }
  if ( tau_DecayMode->at(tauIndex)==11 )
    {
      tau_FR = h_tauFR_11->Eval(tauPt);
      frWeight = tau_FR/(1-tau_FR);
    }
  return frWeight;
}
void mutau_analyzer::fillHist( string histNumber , int muIndex, int tauIndex, float event_weight){
  string hNumber = histNumber;
  plotFill("muPt_"+hNumber,  muPt->at(muIndex) , 30 , 20.0 , 80.0,  event_weight);
  plotFill("muEta_"+hNumber, muEta->at(muIndex), 48, -2.4, 2.4,  event_weight);
  plotFill("muPhi_"+hNumber, muPhi->at(muIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("muDz_"+hNumber,  muDz->at(muIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("muD0_"+hNumber,  muD0->at(muIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("muonID_"+hNumber,muIDbit->at(muIndex)>>1&1, 4, -2, 2,  event_weight); // muonID
  float relMuIso = ( muPFChIso->at(muIndex) + max( muPFNeuIso->at(muIndex) + muPFPhoIso->at(muIndex) - 0.5 *muPFPUIso->at(muIndex) , 0.0 )) / (muPt->at(muIndex));
  plotFill("relMuIso_"+hNumber, relMuIso, 15, 0, 0.3,  event_weight);
  plotFill("muCharge_"+hNumber, muCharge->at(muIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30.0 , 80.0,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 45, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 6, -3, 3,  event_weight);
  //plotFill("tauIso_"+hNumber, tau_IDbits->at(tauIndex)>>16&1, 6, -3, 3,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byVVVLooseDeepTau2017v2p1VSe->at(tauIndex), 4, 0, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byTightDeepTau2017v2p1VSmu->at(tauIndex), 4, 0, 2 ,  event_weight);
  //plotFill("tauAntiEle_"+hNumber, tau_IDbits->at(tauIndex)>>4&1, 8, -2, 2,  event_weight );
  //plotFill("tauAntiMu_"+hNumber,  tau_IDbits->at(tauIndex)>>3&1, 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(muPhi->at(muIndex), muEta->at(muIndex), tau_Phi->at(tauIndex), tau_Eta->at(tauIndex));
  plotFill("deltaR_"+hNumber, deltaR , 40, 0, 6,  event_weight);
  
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(muIndex, tauIndex);
  //if(jetCand.size()>0)
  plotFill("nJet_"+hNumber, jetCand.size() , 6, 0, 6,  event_weight);
  
  plotFill("met_"+hNumber, pfMET , 20, 0, 200,  event_weight);
  
  double mT_muMet = TMass_F((muPt->at(muIndex)),(muPhi->at(muIndex)),pfMET,pfMETPhi  );
  plotFill("mT_muMet_"+hNumber, mT_muMet , 20, 0, 200,  event_weight);

  TLorentzVector myTau; 
  myTau.SetPtEtaPhiE(tau_Pt->at(tauIndex),tau_Eta->at(tauIndex),tau_Phi->at(tauIndex), tau_Energy->at(tauIndex));
  TLorentzVector myMu; 
  myMu.SetPtEtaPhiE(muPt->at(muIndex),muEta->at(muIndex),muPhi->at(muIndex), muE->at(muIndex));
  double visMass_mutau = VisMass_F(myTau, myMu);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  //double HiggsPt = pTvecsum_F(muPt->at(muIndex),tau_Pt->at(tauIndex),muPhi->at(muIndex),tau_Phi->at(tauIndex) );
  double HiggsPt = pTvecsum_F(myMu, myTau, myMet);
  plotFill("higgsPt_"+hNumber,HiggsPt , 40, 0, 400,  event_weight);
  double tot_tr_mass = TotTMass_F(myMu, myTau, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);
  
  int triggerBin=0;
  if( HLTEleMuX>>21&1 == 1 ) triggerBin=1;
  if( HLTEleMuX>>60&1 == 1 ) triggerBin=2;
  if( HLTTau>>0&1 == 1 )     triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 4, 0, 4,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}
void mutau_analyzer::fillHist( string histNumber , TLorentzVector muonP4, TLorentzVector tauP4, int muIndex, int tauIndex, float event_weight){
  string hNumber = histNumber;
  plotFill("muPt_"+hNumber,  muonP4.Pt() , 30 , 20.0 , 80.0,  event_weight);
  plotFill("muEta_"+hNumber, muonP4.Eta(), 48, -2.4, 2.4,  event_weight);
  plotFill("muPhi_"+hNumber, muonP4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("muDz_"+hNumber,  muDz->at(muIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("muD0_"+hNumber,  muD0->at(muIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("muonID_"+hNumber,muIDbit->at(muIndex)>>1&1, 4, -2, 2,  event_weight); // muonID
  float relMuIso = ( muPFChIso->at(muIndex) + max( muPFNeuIso->at(muIndex) + muPFPhoIso->at(muIndex) - 0.5 *muPFPUIso->at(muIndex) , 0.0 )) / (muPt->at(muIndex));
  plotFill("relMuIso_"+hNumber, relMuIso, 15, 0, 0.3,  event_weight);
  plotFill("muCharge_"+hNumber, muCharge->at(muIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 45, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 6, -3, 3,  event_weight);
  //plotFill("tauIso_"+hNumber, tau_IDbits->at(tauIndex)>>16&1, 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byVVVLooseDeepTau2017v2p1VSe->at(tauIndex), 4, 0, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byTightDeepTau2017v2p1VSmu->at(tauIndex), 4, 0, 2 ,  event_weight);
  //plotFill("tauAntiEle_"+hNumber, tau_IDbits->at(tauIndex)>>4&1, 8, -2, 2,  event_weight );
  //plotFill("tauAntiMu_"+hNumber,  tau_IDbits->at(tauIndex)>>3&1, 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(muPhi->at(muIndex), muEta->at(muIndex), tau_Phi->at(tauIndex), tau_Eta->at(tauIndex));
  plotFill("deltaR_"+hNumber, deltaR , 40, 0, 6,  event_weight);
  
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(muIndex, tauIndex);
  //if(jetCand.size()>0)
  plotFill("nJet_"+hNumber, jetCand.size() , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, pfMET , 20, 0, 200,  event_weight);
  
  double mT_muMet = TMass_F((muPt->at(muIndex)),(muPhi->at(muIndex)),pfMET,pfMETPhi  );
  plotFill("mT_muMet_"+hNumber, mT_muMet , 20, 0, 200,  event_weight);

  double visMass_mutau = VisMass_F(muonP4, tauP4);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  //double HiggsPt = pTvecsum_F(muPt->at(muIndex),tau_Pt->at(tauIndex),muPhi->at(muIndex),tau_Phi->at(tauIndex) );
  double HiggsPt = pTvecsum_F(muonP4, tauP4, myMet);
  plotFill("higgsPt_"+hNumber,HiggsPt , 40, 0, 400,  event_weight);
  
  double tot_tr_mass = TotTMass_F(muonP4, tauP4, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);
  
  int triggerBin=0;
  if( HLTEleMuX>>21&1 == 1 ) triggerBin=1;
  if( HLTEleMuX>>60&1 == 1 ) triggerBin=2;
  if( HLTTau>>0&1 == 1 )     triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 4, 0, 4,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
   
}
void mutau_analyzer::fillHist_dyll( string histNumber , int mu1Index, int mu2Index, int tauIndex, float event_weight){
  string hNumber = histNumber;
  plotFill("muPt_"+hNumber,  muPt->at(mu1Index) , 30 , 20 , 80,  event_weight);
  plotFill("muEta_"+hNumber, muEta->at(mu1Index), 48, -2.4, 2.4,  event_weight);
  plotFill("muPhi_"+hNumber, muPhi->at(mu1Index), 30, -3.14, 3.14,  event_weight);
  plotFill("muDz_"+hNumber,  muDz->at(mu1Index), 20, -0.2, 0.2,  event_weight);
  plotFill("muD0_"+hNumber,  muD0->at(mu1Index), 48, -0.06, 0.06,  event_weight);
  plotFill("muonID_"+hNumber,muIDbit->at(mu1Index)>>1&1, 4, -2, 2,  event_weight); // muonID
  float relMuIso = ( muPFChIso->at(mu1Index) + max( muPFNeuIso->at(mu1Index) + muPFPhoIso->at(mu1Index) - 0.5 *muPFPUIso->at(mu1Index) , 0.0 )) / (muPt->at(mu1Index));
  plotFill("relMuIso_"+hNumber, relMuIso, 15, 0, 0.3,  event_weight);
  plotFill("muCharge_"+hNumber, muCharge->at(mu1Index), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 45, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 6, -3, 3,  event_weight);
  //plotFill("tauIso_"+hNumber, tau_IDbits->at(tauIndex)>>16&1, 6, -3, 3,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byVVVLooseDeepTau2017v2p1VSe->at(tauIndex), 4, 0, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byTightDeepTau2017v2p1VSmu->at(tauIndex), 4, 0, 2 ,  event_weight);
  //plotFill("tauAntiEle_"+hNumber, tau_IDbits->at(tauIndex)>>4&1, 8, -2, 2,  event_weight );
  //plotFill("tauAntiMu_"+hNumber,  tau_IDbits->at(tauIndex)>>3&1, 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(muPhi->at(mu1Index), muEta->at(mu1Index), muPhi->at(mu2Index), muEta->at(mu2Index));
  plotFill("deltaR_"+hNumber, deltaR , 40, 0, 6,  event_weight);
  
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(mu1Index, tauIndex);
  //if(jetCand.size()>0)
  plotFill("nJet_"+hNumber, jetCand.size() , 6, 0, 6,  event_weight);
  
  plotFill("met_"+hNumber, pfMET , 20, 0, 200,  event_weight);
  
  double mT_muMet = TMass_F((muPt->at(mu1Index)),(muPhi->at(mu1Index)),pfMET,pfMETPhi  );
  plotFill("mT_muMet_"+hNumber, mT_muMet , 20, 0, 200,  event_weight);

  TLorentzVector myTau; 
  myTau.SetPtEtaPhiE(tau_Pt->at(tauIndex),tau_Eta->at(tauIndex),tau_Phi->at(tauIndex), tau_Energy->at(tauIndex));
  TLorentzVector myMu1; 
  myMu1.SetPtEtaPhiE(muPt->at(mu1Index),muEta->at(mu1Index),muPhi->at(mu1Index), muE->at(mu1Index));
  TLorentzVector myMu2;
  myMu2.SetPtEtaPhiE(muPt->at(mu2Index),muEta->at(mu2Index),muPhi->at(mu2Index), muE->at(mu2Index));
  double visMass_mutau = VisMass_F(myMu1, myMu2);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  //double HiggsPt = pTvecsum_F(muPt->at(mu1Index), muPt->at(mu2Index),muPhi->at(mu1Index),muPhi->at(mu2Index) );
  double HiggsPt = pTvecsum_F(myMu1, myTau, myMet);
  plotFill("higgsPt_"+hNumber,HiggsPt , 40, 0, 400,  event_weight);
  double tot_tr_mass = TotTMass_F(myMu1, myMu2, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);
  
  int triggerBin=0;
  if( HLTEleMuX>>21&1 == 1 ) triggerBin=1;
  if( HLTEleMuX>>60&1 == 1 ) triggerBin=2;
  if( HLTTau>>0&1 == 1 )     triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 4, 0, 4,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}
