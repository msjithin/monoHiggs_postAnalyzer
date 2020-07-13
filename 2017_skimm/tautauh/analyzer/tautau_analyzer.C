/////tautau_analyzer.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom tautau_analyzer analyze
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
#define tautau_analyzer_cxx
#include "tautau_analyzer.h"
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
#include "ComputeFF2018/FFcode/interface/ApplyFF.h"

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
  
  tautau_analyzer t(argv[1],argv[2], isMC);
  t.Loop(maxEvents,reportEvery, SampleName , isMC);
  //delete myMap1;
  cout<<" Outpt written to "<<outputfile<<endl;
  sw.Stop();
  sw.Print();
  return 0;
}

void tautau_analyzer::Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
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
  std::vector<int> tau1Cand;       tau1Cand.clear();
  std::vector<int> tau2Cand;       tau2Cand.clear();
  std::vector<int> jetCand;       jetCand.clear();
  std::vector<int> higgsCand;     higgsCand.clear();
  std::vector<int> reco_tau1;       reco_tau1.clear();
  std::vector<int> reco_tau2;       reco_tau2.clear();
    
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
   TLorentzVector tau1P4;
   TLorentzVector tau2P4;
   TFile frawff("ComputeFF2018/ff_files_tt_2017/uncorrected_fakefactors_tt.root");
   TF1* ff_qcd_0jet=(TF1*) frawff.Get("rawFF_tt_qcd_0jet");
   TF1* ff_qcd_1jet=(TF1*) frawff.Get("rawFF_tt_qcd_1jet");
   TF1* ff_w_0jet=(TF1*) frawff.Get("rawFF_tt_w_0jet");
   TF1* ff_w_1jet=(TF1*) frawff.Get("rawFF_et_w_1jet");
   TF1* ff_tt_0jet=(TF1*) frawff.Get("mc_rawFF_et_tt");

   TFile fmvisclosure ("ComputeFF2018/ff_files_tt_2017/FF_corrections_1.root");
   TF1* mvisclosure_qcd=(TF1*) fmvisclosure.Get("closure_mvis_tt_qcd");
   TF1* mvisclosure_w=(TF1*) fmvisclosure.Get("closure_mvis_tt_w");
   TF1* mvisclosure_tt=(TF1*) fmvisclosure.Get("closure_mvis_tt_ttmc");

   TFile fosssclosure ("ComputeFF2018/ff_files_tt_2017/FF_QCDcorrectionOSSS_tt.root");
   TF1* osssclosure_qcd=(TF1*) fosssclosure.Get("closure_mvis_tt_qcd");
   TF1* mtclosure_w=(TF1*) fosssclosure.Get("closure_mt_tt_w");


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
       tau1Cand.clear();   tau2Cand.clear();
       jetCand.clear();
       higgsCand.clear();
       reco_tau1.clear();   reco_tau2.clear();


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
       bool tauTriggerFilterMatch=false;
       
       bool tauTau_selector=false;
       bool tauDoubleCount=false;
       numberOfEvents+=weight;
       weight=inspected_event_weight;
       if (isMC) genMatching = gen_matching();
       else genMatching = 0;
       //cout<<"genMatch = "<<genMatching<<endl;
       if(debug)cout<<"this worked Line 314"<<endl;
       
       if(isMC) weight=inspected_event_weight;
       else weight=1.0;
       int leading_tau = -1;     float leading_tPt=0;
       int subleading_tau = -1;  float subleading_tPt=0;
       
       if(isMC)
	 pileup_sf = h_pileup->GetBinContent(h_pileup->GetXaxis()->FindBin(puTrue->at(0)));
       weight = weight*pileup_sf;
       if( isGoodVtx==false ) continue;
       //if( noisyJet2017()==true ) continue;
       //if( found_DYjet_sample && !(genMatching<5))
       //	 continue;

       if( found_DYjet_sample && hasGenTau())
	 tauTau_selector=true;
       else if( found_DYjet_sample && !hasGenTau() )
	 tauTau_selector=false;
       else if ( !found_DYjet_sample )
	 tauTau_selector=true;
       
       tauTau_selector=true;
       /////Trigger bit selection
       // if(HLTEleMuX>>21&1 == 1 || HLTEleMuX>>60&1 == 1 )
       // 	 passSingleTriggerPaths=true;
       if( HLTTau>>6&1==1 
	   //|| HLTTau>>10&1==1 
	   || HLTTau>>7&1==1 
	   //|| HLTTau>>11&1==1
	   || HLTTau>>5&1==1
	   //|| HLTTau>>12&1==1
	   )
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
       
       tau1Cand.clear();   tau2Cand.clear();
       event_weight=weight;
       
       ////// signal region -  isolated begin
       if(tauTau_selector)
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
		   tau1Cand = getTauCand(40,2.1);  ///// muons selected 
		   if( tau1Cand.size() >0 ) 
		     { 
		       nGoodMuonPassed+=event_weight;
		       if(debug)cout<<"this worked Line 284"<<endl;
		       tau2Cand = getTau2Cand(40, 2.1, tau1Cand[0]);
		       if( tau2Cand.size()>0 ) 
			 {
			   nGoodTauPassed+=event_weight;
			   if(debug)cout<<"this worked Line 305"<<endl;
			   
			   reco_tau1.clear();reco_tau2.clear();
			   reco_tau1=tau1Cand; reco_tau2=tau2Cand;
			   //cout<<"tau 1 pt: "<<tau_Pt->at(reco_tau1[0])<< "   tau 2 pt: "<<tau_Pt->at(reco_tau2[0])<<endl;
			   if ( MatchTriggerFilter(reco_tau1[0], reco_tau2[0]) )
			     {
			       if ( tau_Charge->at(reco_tau1[0]) * tau_Charge->at(reco_tau2[0]) < 0  ) 
				 {
				   nGoodMuTauPassed+=event_weight;
				   tau1P4.SetPtEtaPhiE(tau_Pt->at(reco_tau1[0]), tau_Eta->at(reco_tau1[0]), tau_Phi->at(reco_tau1[0]), tau_Energy->at(reco_tau1[0]));
				   tau2P4.SetPtEtaPhiE(tau_Pt->at(reco_tau2[0]), tau_Eta->at(reco_tau2[0]), tau_Phi->at(reco_tau2[0]), tau_Energy->at(reco_tau2[0]));
				   
				   if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
				   
				   if(debug)cout<<" sf : "<<getScaleFactors( reco_tau1[0] , reco_tau2[0] , false , isMC , debug ) <<endl;
				   if (isMC) event_weight = event_weight * getScaleFactors( reco_tau1[0] , reco_tau2[0], false , isMC , debug );
				   afterSF4+=event_weight;
				   if( thirdLeptonVeto() < 0 )
				     {
				       nPassedThirdLepVeto+=event_weight;
				       
				       if( passBjetVeto() == true)
					 {
					   nPassedBjetVeto+=event_weight;
					   
					   double deltaR = delta_R(  tau1P4.Phi(), tau1P4.Eta(), tau2P4.Phi(), tau2P4.Eta());
					   if(deltaR > 0.5 )
					     {
					       nDeltaRPassed+=event_weight;
					       if(isMC==false)event_weight=1.0;
					       if(debug)cout<<"this worked Line 374"<<endl;
					       fillHist("5", tau1P4, tau2P4, reco_tau1[0], reco_tau2[0], event_weight);
					       double mT_muMet = TMass_F( tau_Pt->at(reco_tau1[0]),  tau_Phi->at(reco_tau1[0]),pfMET,pfMETPhi  );
					       if( mT_muMet < 50)
						 {
						   fillHist("6", tau1P4, tau2P4, reco_tau1[0], reco_tau2[0], event_weight);
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
	   tau1Cand.clear();   tau2Cand.clear();
	   if(metFilters==0)
	     {
	       
	       if(isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	       nMETFiltersPassed_fr+=event_weight;
	       if(  passSingleTriggerPaths || passCrossTrigger  )
		 {
		   nSingleTrgPassed_fr+=event_weight;
		   if(debug)cout<<"trigger selected"<<endl;
		   tau1Cand = getTauCand(40,2.1);  ///// muons selected 
		   if( tau1Cand.size() >0 ) 
		     { 
		       nGoodMuonPassed_fr+=event_weight;
		       if(debug)cout<<"this worked Line 284"<<endl;
		       tau2Cand = getAISRTau2Cand(40, 2.1, tau1Cand[0]);
		       if( tau2Cand.size()>0 ) 
			 {
			   nGoodTauPassed_fr+=event_weight;
			   if(debug)cout<<"this worked Line 305"<<endl;
			   tauDoubleCount=true;
			   reco_tau1.clear();reco_tau2.clear();
			   reco_tau1=tau1Cand; reco_tau2=tau2Cand;
			   //cout<<"tau 1 pt: "<<tau_Pt->at(reco_tau1[0])<< "   tau 2 pt: "<<tau_Pt->at(reco_tau2[0])<<endl;
			   if ( MatchTriggerFilter(reco_tau1[0], reco_tau2[0]) )
			     {
			       if ( tau_Charge->at(reco_tau1[0]) * tau_Charge->at(reco_tau2[0]) < 0  ) 
				 {
				   nGoodMuTauPassed_fr+=event_weight;
				   tau1P4.SetPtEtaPhiE(tau_Pt->at(reco_tau1[0]), tau_Eta->at(reco_tau1[0]), tau_Phi->at(reco_tau1[0]), tau_Energy->at(reco_tau1[0]));
				   tau2P4.SetPtEtaPhiE(tau_Pt->at(reco_tau2[0]), tau_Eta->at(reco_tau2[0]), tau_Phi->at(reco_tau2[0]), tau_Energy->at(reco_tau2[0]));
				   double mT_muMet = TMass_F( tau1P4.Pt(), tau1P4.Phi(),pfMET,pfMETPhi  );
				   std::vector<int> jetCand;       jetCand.clear();
				   jetCand=getJetCand(reco_tau1[0], reco_tau2[0]);
				   
				   float my_fakefactor = get_ff(tau1P4.Pt(), mT_muMet, VisMass_F(tau1P4, tau2P4), jetCand.size(), 1.0, 0, 0, ff_qcd_0jet, ff_qcd_1jet, 0, 0, 0, mvisclosure_qcd, 0, 0, 0, osssclosure_qcd);
				   cout<<"my_fakefactor = "<<my_fakefactor<<endl;
				   if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
				   
				   event_weight = event_weight* getFR(reco_tau2[0]) * 0.5;
				   //event_weight = event_weight* getFR(reco_tau2[0]);
				   /////
				   if(debug)cout<<" fake bkg sf : "<<getScaleFactors(  reco_tau1[0] , reco_tau2[0] , true , isMC , debug ) <<endl;
				   if(isMC) event_weight = event_weight * getScaleFactors(  reco_tau1[0] , reco_tau2[0] , true , isMC , debug );
				   /////
				   if( thirdLeptonVeto() < 0 )
				     {
				       nPassedThirdLepVeto_fr+=event_weight;
				       
				       if( passBjetVeto() == true)
					 {
					   nPassedBjetVeto_fr+=event_weight;
					   
					   double deltaR = delta_R(  tau1P4.Phi(), tau1P4.Eta(), tau2P4.Phi(), tau2P4.Eta());
					   if(deltaR > 0.5 )
					     {
					       nDeltaRPassed_fr+=event_weight;
					       if(isMC==false)event_weight=1.0;
					       if(debug)cout<<"this worked Line 374"<<endl;
					       fillHist("5_fr", tau1P4, tau2P4, reco_tau1[0], reco_tau2[0], event_weight);
					       double mT_muMet = TMass_F( tau_Pt->at(reco_tau1[0]),  tau_Phi->at(reco_tau1[0]),pfMET,pfMETPhi  );
					       if( mT_muMet < 50)
						 {
						   fillHist("6_fr", tau1P4, tau2P4, reco_tau1[0], reco_tau2[0], event_weight);
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
	   
	   event_weight=weight;
	   tau1Cand.clear();   tau2Cand.clear();
	   
	   if(metFilters==0)
	     {
	       
	       if(isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	       nMETFiltersPassed_fr+=event_weight;
	       if(  passSingleTriggerPaths || passCrossTrigger  )
		 {
		   nSingleTrgPassed_fr+=event_weight;
		   if(debug)cout<<"trigger selected"<<endl;
		   if( tauDoubleCount==false)
		     tau1Cand = getAISRTauCand(40,2.1);  ///// taus selected 
		   else
		     tau1Cand = getTauCand(40,2.1);
		   if( tau1Cand.size() >0 ) 
		     { 
		       nGoodMuonPassed_fr+=event_weight;
		       if(debug)cout<<"this worked Line 284"<<endl;
		       if( tauDoubleCount==false)
			 tau2Cand = getTau2Cand(40, 2.1, tau1Cand[0]);
		       else
			 tau2Cand = getAISRTau2Cand(40, 2.1, tau1Cand[0]);
		       if( tau2Cand.size()>0 ) 
			 {
			   nGoodTauPassed_fr+=event_weight;
			   if(debug)cout<<"this worked Line 305"<<endl;
			   
			   reco_tau1.clear();reco_tau2.clear();
			   reco_tau1=tau1Cand; reco_tau2=tau2Cand;
			   //cout<<"tau 1 pt: "<<tau_Pt->at(reco_tau1[0])<< "   tau 2 pt: "<<tau_Pt->at(reco_tau2[0])<<endl;
			   if ( MatchTriggerFilter(reco_tau1[0], reco_tau2[0]) )
			     {
			       if ( tau_Charge->at(reco_tau1[0]) * tau_Charge->at(reco_tau2[0]) < 0  ) 
				 {
				   nGoodMuTauPassed_fr+=event_weight;
				   tau1P4.SetPtEtaPhiE(tau_Pt->at(reco_tau1[0]), tau_Eta->at(reco_tau1[0]), tau_Phi->at(reco_tau1[0]), tau_Energy->at(reco_tau1[0]));
				   tau2P4.SetPtEtaPhiE(tau_Pt->at(reco_tau2[0]), tau_Eta->at(reco_tau2[0]), tau_Phi->at(reco_tau2[0]), tau_Energy->at(reco_tau2[0]));
				   if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
				   if( tauDoubleCount==false)
				     event_weight = event_weight* getFR(reco_tau1[0])*0.5;
				   else
				     event_weight = event_weight* getFR(reco_tau2[0])*0.5;
				   /////
				   if(debug)cout<<" fake bkg sf : "<<getScaleFactors(  reco_tau1[0] , reco_tau2[0] , true , isMC , debug ) <<endl;
				   if(isMC) event_weight = event_weight * getScaleFactors(  reco_tau1[0] , reco_tau2[0] , true , isMC , debug );
				   /////
				   if( thirdLeptonVeto() < 0 )
				     {
				       nPassedThirdLepVeto_fr+=event_weight;
				       
				       if( passBjetVeto() == true)
					 {
					   nPassedBjetVeto_fr+=event_weight;
					   
					   double deltaR = delta_R(  tau1P4.Phi(), tau1P4.Eta(), tau2P4.Phi(), tau2P4.Eta());
					   if(deltaR > 0.5 )
					     {
					       nDeltaRPassed_fr+=event_weight;
					       if(isMC==false)event_weight=1.0;
					       if(debug)cout<<"this worked Line 374"<<endl;
					       fillHist("5_fr", tau1P4, tau2P4, reco_tau1[0], reco_tau2[0], event_weight);
					       double mT_muMet = TMass_F( tau_Pt->at(reco_tau1[0]),  tau_Phi->at(reco_tau1[0]),pfMET,pfMETPhi  );
					       if( mT_muMet < 50)
						 {
						   fillHist("6_fr", tau1P4, tau2P4, reco_tau1[0], reco_tau2[0], event_weight);
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

void tautau_analyzer::BookHistos(const char* file1, const char* file2)
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


void tautau_analyzer::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{
  
  //  h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);

}




std::vector<int> tautau_analyzer::getTauCand(double tauPtCut, double tauEtaCut){
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
      //if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;

      if(  HLTTau>>6&1 == 1 
      	   || HLTTau>>7&1 == 1
	   || HLTTau>>5&1 == 1
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
std::vector<int> tautau_analyzer::getAISRTauCand(double tauPtCut, double tauEtaCut){
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
      //if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1 )tau_reject=true;
      
      if(  HLTTau>>6&1 == 1 
	   || HLTTau>>7&1 == 1
	   || HLTTau>>5&1 == 1
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
std::vector<int> tautau_analyzer::getTau2Cand(double tauPtCut, double tauEtaCut, int tau1Index){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over taus      
  for(int iTau=0;iTau<nTau;iTau++)
    {
      if( iTau==tau1Index ) continue;  
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
      //if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;

      if(  HLTTau>>6&1 == 1 
	   || HLTTau>>7&1 == 1
	   || HLTTau>>5&1 == 1
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
std::vector<int> tautau_analyzer::getAISRTau2Cand(double tauPtCut, double tauEtaCut, int tau1Index){
  std::vector<int> tmpCand;  tmpCand.clear();
  for(int iTau=0;iTau<nTau;iTau++) //Loop over taus
    {
      if( iTau==tau1Index ) continue;
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
      //if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1 )tau_reject=true;
      if(  HLTTau>>6&1 == 1 
       	   || HLTTau>>7&1 == 1
	   || HLTTau>>5&1 == 1
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
std::vector<int> tautau_analyzer::getJetCand(int tau1Index, int tau2Index){
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
                  
      double dr_jetTau1=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , tau_Phi->at(tau1Index), tau_Eta->at(tau1Index) );
      double dr_jetTau2=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , tau_Phi->at(tau2Index), tau_Eta->at(tau2Index) );
      if( dr_jetTau1>0.5 && dr_jetTau2>0.5 )
	drPassed=true;
	  
      if( (kinematic50 || kinematic30 ) && drPassed==true &&  jetPUFullID->at(iJet)>>1&1==1)
	tmpCand.push_back(iJet);
    }
  return tmpCand;
}
//The noisy jets are defined as: 20 < pt < 50 && abs(eta) > 2.65 && abs(eta) < 3.139. 
bool tautau_analyzer::noisyJet2017(){

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
int tautau_analyzer::thirdLeptonVeto(){
  std::vector<int> tmpCandMu; tmpCandMu.clear();
  std::vector<int> tmpCandEl; tmpCandEl.clear();
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
	tmpCandEl.push_back(iEle);
      }                                                           
    }             
  for(int iMu=0; iMu < nMu;iMu++)
    {
      bool kinematic = false;
      if( (*muPt)[iMu] > 10.0  && fabs((*muEta)[iMu])< 2.4 && (*muD0)[iMu] < 0.045 && (*muDz)[iMu] < 0.2 ) kinematic = true;
      bool muonId =false;
      if( muIDbit->at(iMu)>>1&1==1) muonId =true;
      bool relative_iso = false;
      float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
      if( relMuIso < 0.3 ) relative_iso = true;
      //if( muIDbit->at(iMu)>>6&1==1) relative_iso =true;
      if(muonId==true && kinematic==true && relative_iso==true){
	tmpCandMu.push_back(iMu);
      }                   
    }          

  if(tmpCandEl.size() > 0 || tmpCandMu.size()>0 ){ thirdLepIndex = 1; thirdLepVeto=false;}
  return thirdLepIndex;
  
}


double tautau_analyzer::dR(int mu_index, int tau_index)
{
  double deltaeta = abs(muEta->at(mu_index) - tau_Eta->at(tau_index));
  double muonPhi = muPhi->at(mu_index);
  double tauPhi = tau_Phi->at(tau_index);

  double deltaphi = DeltaPhi(muonPhi, tauPhi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}

double tautau_analyzer::delta_R(float phi1, float eta1, float phi2, float eta2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar   = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}



double tautau_analyzer::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
}

float tautau_analyzer::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
  return  sqrt(2.0*LepPt*met*(1.0-cos(DeltaPhi(LepPhi, metPhi))));
  //return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

float tautau_analyzer::TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float totalTMass = (a + b+ met).M();
  return totalTMass;
}


float tautau_analyzer::VisMass_F(TLorentzVector a, TLorentzVector b){
  float visibleMass = (a + b).M();
  return visibleMass;
}

float tautau_analyzer::pTvecsum_F(float pt1, float pt2, float phi1, float phi2) {
  float pt_vecSum = sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2), 2) + pow(pt1*sin(phi1) + pt2*sin(phi2), 2));
  return pt_vecSum;
}
float tautau_analyzer::pTvecsum_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float pt_vecSum = (a + b+ met).Pt();
  return pt_vecSum;
}

bool tautau_analyzer::passBjetVeto()
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
      return false;
  }
  else if(tmpJetCand.size() >= 2){
    // atleast 2 jets ==> events pass loose
    if( (jetDeepCSVTags_b->at(tmpJetCand[0]) + jetDeepCSVTags_bb->at(tmpJetCand[0])) > 0.1522   )
      return false;
  }
  return true;
}
int tautau_analyzer::gen_matching(){
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
std::vector<int> tautau_analyzer::getGenMu(){
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
bool tautau_analyzer::hasGenTau(){
  bool found_genTau=false;
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>5&1==1) {  found_genTau=true;}
  }
  return found_genTau;
}
float tautau_analyzer::exponential(float x,float a,float b,float c) {
  return a * TMath::Exp(-b * x) + c;
}
double tautau_analyzer::getScaleFactors( int tau1Index , int tau2Index, bool fakeBkg , bool isMC, bool debug)
{
  double rv_sf=1.0;
  double sf_tauidSF_m = 1.0;
  double sf_tauidSF_vvvl = 1.0;
  double sf_tauesSF = 1.0;
  double sf_fakeEle = 1.0; double sf_fakeMu = 1.0;
  double sf_taufesSF = 1.0;
  int genMatchTau = 0;
  if(isMC) genMatchTau = gen_matching();

  sf_tauidSF_m = h_tauidSF_m->GetBinContent(h_tauidSF_m->GetXaxis()->FindBin(tau_DecayMode->at(tau1Index)));
  sf_tauidSF_vvvl = h_tauidSF_vvvl->GetBinContent(h_tauidSF_vvvl->GetXaxis()->FindBin(tau_DecayMode->at(tau1Index)));
  sf_tauesSF = h_tauesSF->GetBinContent(h_tauesSF->GetXaxis()->FindBin(tau_DecayMode->at(tau1Index)));

  //if(genMatchTau==1 || genMatchTau==3)
  //sf_fakeEle = h_tauFakeEleSF->GetBinContent(h_tauFakeEleSF->GetXaxis()->FindBin(abs(tau_Eta->at(tauIndex))));
  //if(genMatchTau==2 || genMatchTau==4)
  //sf_fakeMu = h_tauFakeMuSF->GetBinContent(h_tauFakeMuSF->GetXaxis()->FindBin(abs(tau_Eta->at(tauIndex))));
  if(genMatchTau==2 || genMatchTau==4){
    if(tau_DecayMode->at(tau1Index)==0)
      {
	if(abs(tau_Eta->at(tau1Index)) < 0.4 ) sf_fakeMu=1.14;
	if(abs(tau_Eta->at(tau1Index)) > 0.4 
	   && abs(tau_Eta->at(tau1Index)) < 0.8 ) sf_fakeMu=1.0;
	if(abs(tau_Eta->at(tau1Index)) > 0.8
           && abs(tau_Eta->at(tau1Index)) < 1.2 ) sf_fakeMu=0.87;
	if(abs(tau_Eta->at(tau1Index)) > 1.2
           && abs(tau_Eta->at(tau1Index)) < 1.7 ) sf_fakeMu=0.52;
	if(abs(tau_Eta->at(tau1Index)) > 1.7
           && abs(tau_Eta->at(tau1Index)) < 2.3 ) sf_fakeMu=1.47;
      }
    if(tau_DecayMode->at(tau1Index)==1)
      {
	if(abs(tau_Eta->at(tau1Index)) > 0.0
           && abs(tau_Eta->at(tau1Index)) < 0.4 ) sf_fakeMu=0.69;
      }
  }
  if(genMatchTau==1 || genMatchTau==3){
    if(tau_DecayMode->at(tau1Index)==0)
      {
	if(abs(tau_Eta->at(tau1Index)) < 1.479 ) sf_fakeEle=0.98;
	if(abs(tau_Eta->at(tau1Index)) > 1.479 ) sf_fakeEle=0.80;
      }
    if(tau_DecayMode->at(tau1Index)==1)
      {
	if(abs(tau_Eta->at(tau1Index)) < 1.479 ) sf_fakeEle=1.07;
        if(abs(tau_Eta->at(tau1Index)) > 1.479 ) sf_fakeEle=0.64;
      }
  }
  if(tau_DecayMode->at(tau1Index)==0 && abs(tau_Eta->at(tau1Index))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(1);
  if(tau_DecayMode->at(tau1Index)==0 && abs(tau_Eta->at(tau1Index))>1.4 )  sf_taufesSF = h_taufesSF->Eval(3);
  if(tau_DecayMode->at(tau1Index)==1 && abs(tau_Eta->at(tau1Index))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(5);
  if(tau_DecayMode->at(tau1Index)==1 && abs(tau_Eta->at(tau1Index))>1.4 )  sf_taufesSF = h_taufesSF->Eval(7);
  
  
  //event_weight=event_weight * sf_muID * sf_IsoEff * sf_muTrg * sf_tauidSF_m * sf_tauesSF * sf_fakeEle * (sf_fakeMu);
  // if(fakeBkg)
  //   rv_sf = sf_muID * sf_IsoEff * sf_muTrg * sf_tauidSF_m * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF ;
  // else
  //   rv_sf = sf_muID * sf_IsoEff * sf_muTrg * sf_tauidSF_m * sf_tauidSF_vvvl * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF ;
  if(fakeBkg)
    rv_sf = sf_tauidSF_m * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF;
  else
    rv_sf = sf_tauidSF_m * sf_tauidSF_vvvl * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF;
  
  return rv_sf;

}

bool tautau_analyzer::MatchTriggerFilter(int tau1Index, int tau2Index)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool passFilter = true;
  bool tauTriggerFilterMatch=false;
  int nTauTriggerFilterMatch=0;
  for(int ifilter=0;ifilter<18;ifilter++)
    {
      if(tauFiredTrgs->at(tau1Index)>>ifilter&1==1)
	{
	  tauTriggerFilterMatch=true;
	  nTauTriggerFilterMatch++;
	}
    }

  if(  HLTTau>>6&1 == 1 
       || HLTTau>>7&1 == 1
       || HLTTau>>5&1 == 1
       )  passFilter=true;
    
  return passFilter;
}



double  tautau_analyzer::getFR(int tauIndex){
  double frWeight=1.0;
  double tau_FR = 1.0;
  double tauPt=0.0;
  if( tau_Pt->at(tauIndex) < 120 )
    tauPt=tau_Pt->at(tauIndex);
  else
    tauPt=119.0;
  // if ( tau_DecayMode->at(tauIndex)==0 )
  //   {
  //     tau_FR = h_tauFR_0->Eval(tauPt);
  //     frWeight = tau_FR/(1-tau_FR);
  //   }
  
  // if ( tau_DecayMode->at(tauIndex)==1 )
  //   {
  //     tau_FR = h_tauFR_1->Eval(tauPt);
  //     frWeight = tau_FR/(1-tau_FR);
  //   }
  
  // if ( tau_DecayMode->at(tauIndex)==10 )
  //   {
  //     tau_FR = h_tauFR_10->Eval(tauPt);
  //     frWeight = tau_FR/(1-tau_FR);
  //   }
  // if ( tau_DecayMode->at(tauIndex)==11 )
  //   {
  //     tau_FR = h_tauFR_11->Eval(tauPt);
  //     frWeight = tau_FR/(1-tau_FR);
  //   }
  if ( tau_DecayMode->at(tauIndex)==0 )
    {
      tau_FR = h_tauFR_0->GetBinContent(h_tauFR_0->GetXaxis()->FindBin(tauPt));
      frWeight = tau_FR/(1-tau_FR);
    }
  
  if ( tau_DecayMode->at(tauIndex)==1 )
    {
      tau_FR = h_tauFR_1->GetBinContent(h_tauFR_1->GetXaxis()->FindBin(tauPt));
      frWeight = tau_FR/(1-tau_FR);
    }
  
  if ( tau_DecayMode->at(tauIndex)==10 )
    {
      tau_FR = h_tauFR_10->GetBinContent(h_tauFR_10->GetXaxis()->FindBin(tauPt));
      frWeight = tau_FR/(1-tau_FR);
    }
  return frWeight;
}
void tautau_analyzer::fillHist( string histNumber , int tau1Index, int tau2Index, float event_weight){
  string hNumber = histNumber;
  plotFill("tau1Pt_"+hNumber,  tau_Pt->at(tau1Index) , 25 , 30 , 80,  event_weight);
  plotFill("tau1Eta_"+hNumber, tau_Eta->at(tau1Index), 45, -2.5, 2.5,  event_weight);
  plotFill("tau1Phi_"+hNumber, tau_Phi->at(tau1Index), 30, -3.14, 3.14,  event_weight);
  plotFill("tau1Iso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tau1Index), 6, -3, 3,  event_weight);
  plotFill("tau1DecayMode_"+hNumber, tau_DecayMode->at(tau1Index) , 12, 0, 12,  event_weight);
  plotFill("tau1Charge_"+hNumber, tau_Charge->at(tau1Index), 8, -2, 2 ,  event_weight);
  plotFill("tau1AntiEle_"+hNumber, tau_byVVVLooseDeepTau2017v2p1VSe->at(tau1Index), 4, 0, 2,  event_weight );
  plotFill("tau1AntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tau1Index), 4, 0, 2 ,  event_weight);

  plotFill("tau2Pt_"+hNumber,  tau_Pt->at(tau2Index) , 25 , 30 , 80,  event_weight);
  plotFill("tau2Eta_"+hNumber, tau_Eta->at(tau2Index), 45, -2.5, 2.5,  event_weight);
  plotFill("tau2Phi_"+hNumber, tau_Phi->at(tau2Index), 30, -3.14, 3.14,  event_weight);
  plotFill("tau2Iso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tau2Index), 6, -3, 3,  event_weight);
  plotFill("tau2DecayMode_"+hNumber, tau_DecayMode->at(tau2Index) , 12, 0, 12,  event_weight);
  plotFill("tau2Charge_"+hNumber, tau_Charge->at(tau2Index), 8, -2, 2 ,  event_weight);
  plotFill("tau2AntiEle_"+hNumber, tau_byVVVLooseDeepTau2017v2p1VSe->at(tau2Index), 4, 0, 2,  event_weight );
  plotFill("tau2AntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tau2Index), 4, 0, 2 ,  event_weight);
  double deltaR = delta_R(tau_Phi->at(tau1Index), tau_Eta->at(tau1Index), tau_Phi->at(tau2Index), tau_Eta->at(tau2Index));
  plotFill("deltaR_"+hNumber, deltaR , 40, 0, 6,  event_weight);
  
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(tau1Index, tau2Index);
  plotFill("nJet_"+hNumber, jetCand.size() , 6, 0, 6,  event_weight);
  
  plotFill("met_"+hNumber, pfMET , 20, 0, 200,  event_weight);
  
  double mT_muMet = TMass_F( tau_Pt->at(tau1Index),  tau_Phi->at(tau1Index),pfMET,pfMETPhi  );
  plotFill("mT_muMet_"+hNumber, mT_muMet , 20, 0, 200,  event_weight);

  TLorentzVector myTau1; 
  myTau1.SetPtEtaPhiE(tau_Pt->at(tau1Index),tau_Eta->at(tau1Index),tau_Phi->at(tau1Index), tau_Energy->at(tau1Index));
  TLorentzVector myTau2; 
  myTau2.SetPtEtaPhiE(tau_Pt->at(tau2Index),tau_Eta->at(tau2Index),tau_Phi->at(tau2Index), tau_Energy->at(tau2Index));

  double visMass_mutau = VisMass_F(myTau1, myTau2);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  //double HiggsPt = pTvecsum_F(muPt->at(muIndex),tau_Pt->at(tauIndex),muPhi->at(muIndex),tau_Phi->at(tauIndex) );
  double HiggsPt = pTvecsum_F(myTau1, myTau2, myMet);
  plotFill("higgsPt_"+hNumber,HiggsPt , 40, 0, 400,  event_weight);
  double tot_tr_mass = TotTMass_F(myTau1, myTau1, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);
  
  int triggerBin=0;
  if(  HLTTau>>6&1 == 1 ||  HLTTau>>10&1 == 1 ) triggerBin=1;
  if(  HLTTau>>7&1 == 1 ||  HLTTau>>11&1 == 1 ) triggerBin=2;
  if( HLTTau>>5&1 == 1 )     triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 4, 0, 4,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}
void tautau_analyzer::fillHist( string histNumber , TLorentzVector tau1P4, TLorentzVector tau2P4, int tau1Index, int tau2Index, float event_weight){
  string hNumber = histNumber;
  
  plotFill("tau1Pt_"+hNumber,  tau_Pt->at(tau1Index) , 25 , 30 , 80,  event_weight);
  plotFill("tau1Eta_"+hNumber, tau_Eta->at(tau1Index), 45, -2.5, 2.5,  event_weight);
  plotFill("tau1Phi_"+hNumber, tau_Phi->at(tau1Index), 30, -3.14, 3.14,  event_weight);
  plotFill("tau1Iso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tau1Index), 6, -3, 3,  event_weight);
  plotFill("tau1DecayMode_"+hNumber, tau_DecayMode->at(tau1Index) , 12, 0, 12,  event_weight);
  plotFill("tau1Charge_"+hNumber, tau_Charge->at(tau1Index), 8, -2, 2 ,  event_weight);
  plotFill("tau1AntiEle_"+hNumber, tau_byVVVLooseDeepTau2017v2p1VSe->at(tau1Index), 4, 0, 2,  event_weight );
  plotFill("tau1AntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tau1Index), 4, 0, 2 ,  event_weight);

  plotFill("tau2Pt_"+hNumber,  tau_Pt->at(tau2Index) , 25 , 30 , 80,  event_weight);
  plotFill("tau2Eta_"+hNumber, tau_Eta->at(tau2Index), 45, -2.5, 2.5,  event_weight);
  plotFill("tau2Phi_"+hNumber, tau_Phi->at(tau2Index), 30, -3.14, 3.14,  event_weight);
  plotFill("tau2Iso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tau2Index), 6, -3, 3,  event_weight);
  plotFill("tau2DecayMode_"+hNumber, tau_DecayMode->at(tau2Index) , 12, 0, 12,  event_weight);
  plotFill("tau2Charge_"+hNumber, tau_Charge->at(tau2Index), 8, -2, 2 ,  event_weight);
  plotFill("tau2AntiEle_"+hNumber, tau_byVVVLooseDeepTau2017v2p1VSe->at(tau2Index), 4, 0, 2,  event_weight );
  plotFill("tau2AntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tau2Index), 4, 0, 2 ,  event_weight);
  double deltaR = delta_R(tau_Phi->at(tau1Index), tau_Eta->at(tau1Index), tau_Phi->at(tau2Index), tau_Eta->at(tau2Index));
  plotFill("deltaR_"+hNumber, deltaR , 40, 0, 6,  event_weight);
  
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(tau1Index, tau2Index);
  //if(jetCand.size()>0)
  plotFill("nJet_"+hNumber, jetCand.size() , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, pfMET , 20, 0, 200,  event_weight);
  
  double mT_muMet = TMass_F( tau_Pt->at(tau1Index),  tau_Phi->at(tau1Index),pfMET,pfMETPhi  );
  plotFill("mT_muMet_"+hNumber, mT_muMet , 20, 0, 200,  event_weight);

  double visMass_mutau = VisMass_F(tau1P4, tau2P4);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  //double HiggsPt = pTvecsum_F(muPt->at(muIndex),tau_Pt->at(tauIndex),muPhi->at(muIndex),tau_Phi->at(tauIndex) );
  double HiggsPt = pTvecsum_F(tau1P4, tau2P4, myMet);
  plotFill("higgsPt_"+hNumber,HiggsPt , 40, 0, 400,  event_weight);
  
  double tot_tr_mass = TotTMass_F(tau1P4, tau2P4, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);
  
  int triggerBin=0;
  if(  HLTTau>>6&1 == 1 ||  HLTTau>>10&1 == 1 ) triggerBin=1;
  if(  HLTTau>>7&1 == 1 ||  HLTTau>>11&1 == 1 ) triggerBin=2;
  if( HLTTau>>5&1 == 1 )     triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 4, 0, 4,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
   
}
