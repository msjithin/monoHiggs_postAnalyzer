/////etau_analyzer.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom etau_analyzer analyze
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
#define etau_analyzer_cxx
#include "etau_analyzer.h"
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
#include "ApplyScaleFactors.h"
#include "commonFunctions.h"
#include "fillEtauHistograms.h"
//#include "TauTriggerSFs2017.cc"
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
  
  etau_analyzer t(argv[1],argv[2], isMC);
  t.Loop(maxEvents,reportEvery, SampleName , isMC);
  //delete myMap1;
  cout<<" Outpt written to "<<outputfile<<endl;
  sw.Stop();
  sw.Print();
  return 0;
}

void etau_analyzer::Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
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
  if(debug) cout<<"******** debugging is on ******************"<<endl;
  double netWeight = 1.0;
  double afterSF1=0;
  double afterSF2=0;     
  double afterSF3=0;     
  double afterSF4=0;     

  if (fChain == 0) return;
  int genMatching=0; 
  int thirdLeptonIndex=-1;
  std::vector<int> muonGen;
  std::vector<int> eleCand;        eleCand.clear();
  std::vector<int> ele2Cand;       ele2Cand.clear();
  std::vector<int> tauCand;        tauCand.clear();
  std::vector<int> reco_ele;      reco_ele.clear(); 
  std::vector<int> reco_tau;      reco_tau.clear(); 
  std::vector<int> reco_ele2;      reco_ele2.clear();
  std::vector<int> jetCand;       jetCand.clear();
  
  TString sample = TString(SampleName);
  int nHiggs = 0;
  bool fill_hist = false;
  bool isMC = false;
  if( _isMC_=="MC" ) { isMC=true; fill_hist=true; }
  else if ( _isMC_=="DATA" ) { isMC=false; fill_hist=false; }
  
  Double_t  Pt_Bins[26]={0.0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
  Double_t  Pt_Bins_highPt[21]={100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
  
  TH1F* h_cutflow=new TH1F("cutflow", "cutflow", 10, 0, 10); h_cutflow->Sumw2();
  TH1F* h_cutflow_n=new TH1F("cutflow_n", "cutflow_n", 8, 0, 8);h_cutflow_n->Sumw2();
  TH1F* h_cutflow_n_fr=new TH1F("cutflow_n_fr", "cutflow_n_fr", 8, 0, 8);h_cutflow_n_fr->Sumw2();
  TH1F* h_cutflow_n_dyll=new TH1F("cutflow_n_dyll", "cutflow_n_dyll", 8, 0, 8);h_cutflow_n_dyll->Sumw2();
  TH1F* h_cutflow_n_dyll_fr=new TH1F("cutflow_n_dyll_fr", "cutflow_n_dyll_fr", 8, 0, 8);h_cutflow_n_dyll_fr->Sumw2();
  //TH1F* h_cutflow_Htt=new TH1F("cutflow_Htt", "cutflow_Htt", 11, 0, 11); h_cutflow_Htt->Sumw2();
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
   TFile* f_nnlo_qcd = TFile::Open("sf_files/2017/RootFiles/theory/lindert_qcd_nnlo_sf.root");
   TFile* f_nlo_qcd  = TFile::Open("sf_files/2017/RootFiles/theory/2017_gen_v_pt_qcd_sf.root");
   TFile* f_qcd_ewk;
   if ( found_Wjet_sample ) {
     f_qcd_ewk = TFile::Open("sf_files/2017/RootFiles/theory/merged_kfactors_wjets.root");
     NLO_QCD_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_qcd_ewk");
     NLO_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_ewk");
     NLO_QCD = (TH1F*)f_nlo_qcd->Get("wjet_dress_monojet");
     NNLO_QCD = (TH1F*)f_nnlo_qcd->Get("evj");
     
   } else if ( found_DYjet_sample ) {
     f_qcd_ewk = TFile::Open("sf_files/2017/RootFiles/theory/merged_kfactors_zjets.root");
     NLO_QCD_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_qcd_ewk");
     NLO_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_ewk");
     f_nlo_qcd = TFile::Open("sf_files/2017/RootFiles/theory/kfac_dy_filter.root");
     NLO_QCD = (TH1F*)f_nlo_qcd->Get("kfac_dy_filter");
     NNLO_QCD = (TH1F*)f_nnlo_qcd->Get("eej");
     
   }

  
   TLorentzVector electronP4;
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
       eleCand.clear();
       tauCand.clear();
       reco_ele.clear(); reco_ele2.clear();
       reco_tau.clear();  
       jetCand.clear();
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
       //double muRC_sf = 1.0; double randomN = gRandom->Rndm();

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
       bool eleTriggerFilterMatch=false;
       bool tauTriggerFilterMatch=false;
       bool eleTau_selector=false;
       bool eleEle_selector=false;

       numberOfEvents+=weight;
       if(isMC) weight=inspected_event_weight;
       else weight=1.0;
       if(isMC)
       	 pileup_sf = h_pileup->GetBinContent(h_pileup->GetXaxis()->FindBin(puTrue->at(0)));
       weight = weight*pileup_sf;
       // if(isMC)
       // 	 weight=weight*prefiringweight;
       if( isGoodVtx==false ) continue;
       //if( noisyJet2017()==true ) continue;
       if( found_DYjet_sample && hasGenTau())
	 eleTau_selector=true;
       else if( found_DYjet_sample && !hasGenTau() )
	 eleTau_selector=false;
       else if ( !found_DYjet_sample )
	 eleTau_selector=true;

       if( found_DYjet_sample )
	 {
	   if(!found_GenMatch(5) && !found_GenMatch(6) )
	     eleEle_selector=true;
	   else
	     eleEle_selector=false;
	 }
       else
	 eleEle_selector=false;

       eleTau_selector=true;
       eleEle_selector=true;
       ////Trigger bit selection
       if( (HLTEleMuX>>3&1 == 1 )
	   || (HLTEleMuX>>61&1 == 1 )
	   || (HLTEleMuX>>5&1 == 1)
	   )
       	 passSingleTriggerPaths=true;
       
       if( ( HLTTau>>1&1 == 1 ) )
       	 passCrossTrigger=true;
       
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
       //weight=weight*kfactor;
       // if(kfactor!=1)cout<<" kfactor = "<<kfactor<<endl;
       // cout<<" weight = "<<weight<<endl;
       if(debug)cout<<"this worked Line 332"<<endl;
              
       /////
       if(debug)cout<<"entry # : "<<jentry<<endl;
       event_weight=weight;
       if(debug)cout<<"reco selections begin"<<endl;
       eleCand.clear(); ele2Cand.clear();  tauCand.clear();
       

       ////// reco selection begin
       if(debug)cout<<"signal region DY->ll -  isolated begin"<<endl;
       ////// DY Z-> ll signal region -  isolated begin
       if(eleEle_selector)
	 {
	   if(metFilters==0 )
	     {
	       if(debug)cout<<"metfilters selected"<<endl;
	       if(isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	       nMETFiltersPassed_dyll+=event_weight;
	       makeTestPlot("a_dyll", 0,0,0,event_weight);
	       if(debug)cout<<"genweight applied"<<endl;
	       if(   passSingleTriggerPaths || passCrossTrigger )
		 {
		   nSingleTrgPassed_dyll+=event_weight;
		   if(debug)cout<<"trigger selected"<<endl;
		   makeTestPlot("b_dyll", 0,0,0,event_weight);
		   eleCand = getEleCand(24,2.1);  ///// ele selected
		   if( eleCand.size() >0 ) 
		     { 
		       nGoodMuonPassed_dyll+=event_weight;
		       if(debug)cout<<"this worked Line 443"<<endl;
		       makeTestPlot("c_dyll", 0,0,0,event_weight);
		       std::vector<int> iElePlus; iElePlus.clear(); 
		       std::vector<int> iEleMinus; iEleMinus.clear();
		       for(int i=0; i<eleCand.size(); i++){
			 if(eleCharge->at(eleCand[i]) < 0) iEleMinus.push_back(eleCand[i]);
			 if(eleCharge->at(eleCand[i]) > 0) iElePlus.push_back(eleCand[i]);
		       }
		       tauCand = getTauCand(30,2.3);
		       if( tauCand.size()>0 
			   //&&  tau_byTightDeepTau2017v2p1VSe->at(tauCand[0])==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(tauCand[0])==1 
			   ) 
			 {
			   nGoodTauPassed_dyll+=event_weight;
			   if(debug)cout<<"this worked Line 424"<<endl;
			   makeTestPlot("d_dyll", 0,0,0,event_weight);
			   if(iElePlus.size()>0 && iEleMinus.size()>0)
			     {
			       if( eleCharge->at(iElePlus[0])*eleCharge->at(iEleMinus[0]) <0 )
				 {
				   reco_ele.clear();reco_tau.clear(); reco_ele2.clear();
				   if( elePt->at(iEleMinus[0]) > elePt->at(iElePlus[0])   ) { reco_ele=iEleMinus; reco_ele2=iElePlus; }
				   else { reco_ele=iElePlus; reco_ele2=iEleMinus;}
				   reco_tau=tauCand;
				   //reco_tau=reco_ele2;
				   makeTestPlot("e_dyll", 0,0,0,event_weight);
				   if ( MatchTriggerFilter(reco_ele[0], reco_tau[0]) )
				     {
				   
				       //if ( eleCharge->at(reco_ele[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
				       {
					 nGoodMuTauPassed_dyll+=event_weight;
					 //electronP4.SetPtEtaPhiE(elePt->at(reco_ele[0]), eleEta->at(reco_ele[0]), elePhi->at(reco_ele[0]), eleE->at(reco_ele[0]));
					 //tauP4.SetPtEtaPhiE(tau_Pt->at(reco_tau[0]), tau_Eta->at(reco_tau[0]), tau_Phi->at(reco_tau[0]), tau_Energy->at(reco_tau[0]));
					 if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
					 
					 if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , false , isMC , debug ) <<endl;
					 if (isMC) event_weight = event_weight * getScaleFactors( reco_ele[0] , reco_tau[0] , false , isMC , debug );
					 afterSF4+=event_weight;
					 makeTestPlot("f_dyll", 0,0,0,event_weight);
					 if( thirdLeptonVeto(reco_ele[0] , reco_tau[0]) < 0 )
					   {
					     nPassedThirdLepVeto_dyll+=event_weight;
					     makeTestPlot("g_dyll", 0,0,0,event_weight);
					     if( passBjetVeto() == true)
					       {
						 nPassedBjetVeto_dyll+=event_weight;
						 makeTestPlot("h_dyll", 0,0,0,event_weight);
						 double deltaR = delta_R(elePhi->at(reco_ele[0]),eleEta->at(reco_ele[0]), elePhi->at(reco_ele2[0]),eleEta->at(reco_ele2[0]));
						 if(deltaR > 0.5 )
						   {
						     nDeltaRPassed_dyll+=event_weight;
						     if(isMC==false)event_weight=1.0;
						     makeTestPlot("i_dyll", 0,0,0,event_weight);
						     if(debug)cout<<"this worked Line 374"<<endl;
						     fillHist_dyll("5_dyll",  reco_ele[0], reco_ele2[0], reco_tau[0], event_weight, isMC);
						     double mT_muMet = TMass_F((elePt->at(reco_ele[0])),(elePhi->at(reco_ele[0])),pfMET,pfMETPhi  );
						     if( mT_muMet < 50 )
						       {
							 fillHist_dyll("6_dyll", reco_ele[0], reco_ele2[0], reco_tau[0], event_weight, isMC);
							 makeTestPlot("j_dyll", 0,0,0,event_weight);
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
	 }
       if(debug)cout<<"signal region -  isolated begin L523"<<endl;       
       
       bool Ztt_selector=false;
       
       ////// signal region -  isolated begin
       event_weight=weight;
       eleCand.clear(); ele2Cand.clear();  tauCand.clear();
       if(eleTau_selector)
	 {
	   if(metFilters==0)
	     {
	       if(debug)cout<<"metfilters selected"<<endl;
	       if (isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	       nMETFiltersPassed+=event_weight;
	       makeTestPlot("a", 0,0,0,event_weight);
	       if(debug)cout<<"genweight applied"<<endl;
	       if( passSingleTriggerPaths || passCrossTrigger  )
		 {
		   nSingleTrgPassed+=event_weight;
		   if(debug)cout<<"trigger selected"<<endl;
		   makeTestPlot("b", 0,0,0,event_weight);
		   eleCand = getEleCand(24,2.1);  ///// ele selected 
		   if( eleCand.size() >0 ) 
		     { 
		       nGoodMuonPassed+=event_weight;
		       if(debug)cout<<"this worked Line 526"<<endl;
		       
		       makeTestPlot("c", 0,0,0,event_weight);
		       tauCand = getTauCand(30,2.3);
		       if( tauCand.size() >0 )
			 {
			   nGoodTauPassed+=event_weight;
			   reco_ele.clear();reco_tau.clear();
			   reco_ele=eleCand; reco_tau=tauCand;
			   makeTestPlot("d", 0,0,0,event_weight);
			   
			   if(found_DYjet_sample){
			     if( passDiElectronVeto(eleCand[0])==true 
				 && !(eVetoZTTp001dxyz(0.001)>1)
				 && !(mVetoZTTp001dxyz(0.001)>0)
				 ) Ztt_selector=true;
			     else Ztt_selector=false;
			   }
			   else
			     Ztt_selector=true;
			   
			   //if(Ztt_selector) 
			     {
			   if (  eleCharge->at(reco_ele[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
			     {
			       nGoodMuTauPassed+=event_weight;
			       if(debug)cout<<"this worked Line 538"<<endl;
			       afterSF1+=event_weight;
			       makeTestPlot("e", 0,0,0,event_weight);
			       if ( MatchTriggerFilter(reco_ele[0], reco_tau[0]) )
				 {
				   if(debug)cout<<"this worked Line 534"<<endl;
				   
				   if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , false , isMC , debug ) <<endl;
				   if (isMC) event_weight = event_weight * getScaleFactors( reco_ele[0] , reco_tau[0] , false , isMC , debug );
				   
				   afterSF4+=event_weight;
				   makeTestPlot("f", 0,0,0,event_weight);
				   if( thirdLeptonVeto(reco_ele[0] , reco_tau[0]) < 0 )
				     {
				       nPassedThirdLepVeto+=event_weight;
				       makeTestPlot("g", 0,0,0,event_weight);
				       if( passBjetVeto() == true)
					 {
					   nPassedBjetVeto+=event_weight;
					   makeTestPlot("h", 0,0,0,event_weight);
					   double deltaR = delta_R(elePhi->at(reco_ele[0]),eleEta->at(reco_ele[0]), tau_Phi->at(reco_tau[0]),  tau_Eta->at(reco_tau[0]));
					   if(deltaR > 0.5 )
					     {
					       nDeltaRPassed+=event_weight;
					       if(isMC==false)event_weight=1.0;
					       if(debug)cout<<"this worked Line 558"<<endl;
					       fillHist("5", reco_ele[0], reco_tau[0], event_weight, isMC);
					       makeTestPlot("i", 0,0,0,event_weight);
					       double mT_eleMet = TMass_F((elePt->at(reco_ele[0])),(elePhi->at(reco_ele[0])),pfMET,pfMETPhi  );
					       if( mT_eleMet < 50 )
						 {
						   fillHist("6", reco_ele[0], reco_tau[0], event_weight, isMC);
						   makeTestPlot("j", 0,0,0,event_weight);
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
       //////// signal region end
       if(debug)cout<<"fake background region - antiisolated begin 625"<<endl;
       ///// fake background region - antiisolated begin
       event_weight=weight;
       eleCand.clear(); tauCand.clear();
       if(metFilters==0)
	 {
	   if (isMC)fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed_fr+=event_weight;
	   makeTestPlot("a_fr", 0,0,0,event_weight);
	   if(  passSingleTriggerPaths || passCrossTrigger  )
	     {
	       nSingleTrgPassed_fr+=event_weight;
	       if(debug)cout<<"trigger selected line 636"<<endl;
	       makeTestPlot("b_fr", 0,0,0,event_weight);
	       eleCand = getEleCand(24,2.1);  ///// ele selected 
	       if( eleCand.size() >0 ) 
		 { 
		   nGoodMuonPassed_fr+=event_weight;
		   makeTestPlot("c_fr", 0,0,0,event_weight);
		   if(debug)cout<<"this worked Line 641"<<endl;
		   tauCand = getAISRTauCand(30,2.3);
		   if( tauCand.size()>0 ) 
		     {
		       nGoodTauPassed_fr+=event_weight;
		       makeTestPlot("d_fr", 0,0,0,event_weight);
		       reco_ele.clear();reco_tau.clear();
		       reco_ele=eleCand; reco_tau=tauCand;
		       if (  eleCharge->at(reco_ele[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
			 {
			   nGoodMuTauPassed_fr+=event_weight;
			   makeTestPlot("e_fr", 0,0,0,event_weight);
			   if ( MatchTriggerFilter(reco_ele[0], reco_tau[0]) )
			     {
			       
			       if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , true , isMC , debug ) <<endl;
			       if (isMC) event_weight = event_weight * getScaleFactors( reco_ele[0] , reco_tau[0] , true , isMC , debug );
			       
			       if (debug==true ) std::cout<<"event_weight =  "<< event_weight<<" event number = "<<jentry <<std::endl;
			       event_weight = event_weight* getFR(reco_tau[0]);
			       makeTestPlot("f_fr", 0,0,0,event_weight);
			       if( thirdLeptonVeto(reco_ele[0] , reco_tau[0]) < 0 )
				 {
				   nPassedThirdLepVeto_fr+=event_weight;
				   makeTestPlot("g_fr", 0,0,0,event_weight);
				   if( passBjetVeto() == true)
				     {
				       nPassedBjetVeto_fr+=event_weight;
				       makeTestPlot("h_fr", 0,0,0,event_weight);
				       double deltaR = delta_R(elePhi->at(reco_ele[0]),eleEta->at(reco_ele[0]), tau_Phi->at(reco_tau[0]),  tau_Eta->at(reco_tau[0]));
				       if(deltaR > 0.5 )
					 {
					   nDeltaRPassed_fr+=event_weight;
					   makeTestPlot("i_fr", 0,0,0,event_weight);
					   if(debug)cout<<"this worked Line 442"<<endl;
					   fillHist("5_fr", reco_ele[0], reco_tau[0], event_weight, isMC);
					   double mT_eleMet = TMass_F((elePt->at(reco_ele[0])),(elePhi->at(reco_ele[0])),pfMET,pfMETPhi  );
					   if( mT_eleMet < 50 )
					     {
					       fillHist("6_fr", reco_ele[0], reco_tau[0], event_weight, isMC);
					       makeTestPlot("j_fr", 0,0,0,event_weight);
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

   // std::cout<<std::setw(20) <<std::right <<"after sf 1 "<<afterSF1<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nGoodMuTauPassed)*100/nGoodTauPassed<<std::endl;
   // std::cout<<std::setw(20) <<std::right <<"after sf 2 "<<afterSF2<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nGoodMuTauPassed)*100/nGoodTauPassed<<std::endl;
   // std::cout<<std::setw(20) <<std::right <<"after sf 3 "<<afterSF3<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nGoodMuTauPassed)*100/nGoodTauPassed<<std::endl;
   // std::cout<<std::setw(20) <<std::right <<"after sf 4 "<<afterSF4<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nGoodMuTauPassed)*100/nGoodTauPassed<<std::endl;

   
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




//---------------------------------------------------                                                                                                                                
// get a electron candiate based on pt eta and isolation                                                                                                                               
//----------------------------------------------------                                                                                                                               

std::vector<int> etau_analyzer::getEleCand(double elePtCut, double eleEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over electrons                                                                     
  for(int iEle=0;iEle<nEle;iEle++)
    {
      bool kinematic = false;
      if( elePt->at(iEle) > elePtCut  
	  && fabs(eleEta->at(iEle))< eleEtaCut 
	  && fabs(eleD0->at(iEle)) < 0.045
	  && fabs(eleDz->at(iEle)) < 0.2
	  && eleMissHits->at(iEle) <= 1 && eleConvVeto->at(iEle)==1
	  ) kinematic = true;
      bool electronId =false;
      if( eleIDbit->at(iEle)>>8&1==1) electronId =true;
      bool relative_iso = false;    
      float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
      if( relEleIso < 0.15 ) relative_iso = true;
      bool trigger = false;
      if( ( HLTEleMuX>>3&1 == 1 && elePt->at(iEle) > 28.0  ) 
	  || ( HLTEleMuX>>61&1 == 1 && elePt->at(iEle) > 33.0  )
	  || ( HLTEleMuX>>5&1 == 1 && elePt->at(iEle) > 36.0 )
	  || ( HLTTau>>1&1 == 1 && elePt->at(iEle) > 25.0  && elePt->at(iEle) < 28.0 && fabs(eleEta->at(iEle))< 2.1 )

	   ) trigger = true;
      if( kinematic && electronId && relative_iso && trigger ){
	tmpCand.push_back(iEle);
      }	
    }                           
  return tmpCand;
  
}

std::vector<int> etau_analyzer::getTauCand(double tauPtCut, double tauEtaCut){
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
      bool trigger = false;
      if( tau_Pt->at(iTau) > tauPtCut 
	  && fabs( tau_Eta->at(iTau))< tauEtaCut 
	  && tau_LeadChargedHadron_dz->at(iTau) < 0.2
	  //&& fabs(tau_Charge->at(iTau))==1
	  )kinematic = true;
      //if( tau_IDbits->at(iTau)>>16&1==1 ) tauIsolation=true;
      if( tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1 ) tauIsolation=true; 
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byTightDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      if( ( HLTEleMuX>>3&1 == 1 ) 
	  || ( HLTEleMuX>>61&1 == 1 )
	  || ( HLTEleMuX>>5&1 == 1 )
	  || ( HLTTau>>1&1 == 1 && tau_Pt->at(iTau) >35 && fabs(tau_Eta->at(iTau)) < 2.1 )
	  
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
std::vector<int> etau_analyzer::getAISRTauCand(double tauPtCut, double tauEtaCut){
  std::vector<int> tmpCand;  tmpCand.clear();
  for(int iTau=0;iTau<nTau;iTau++) //Loop over taus
    {
      bool kinematic = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool mutau_separation=false;
      bool newDecayModeFinding=false;
      bool tau_reject=false;
      bool trigger = false;
      if( tau_Pt->at(iTau) > tauPtCut 
	  && fabs( tau_Eta->at(iTau))< tauEtaCut 
	  && tau_LeadChargedHadron_dz->at(iTau) < 0.2
	  && fabs(tau_Charge->at(iTau))==1
  	  )kinematic = true;
      if(  tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 && tau_byMediumDeepTau2017v2p1VSjet->at(iTau)!=1 ) tauIsolation=true;
      //if( tau_IDbits->at(iTau)>>13&1==1 && !(tau_IDbits->at(iTau)>>16&1==1) ) tauIsolation=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byTightDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      if( ( HLTEleMuX>>3&1 == 1 )
	  || ( HLTEleMuX>>61&1 == 1 )
          || ( HLTEleMuX>>5&1 == 1 )
          || ( HLTTau>>1&1 == 1 && tau_Pt->at(iTau) >35 && abs(tau_Eta->at(iTau)) < 2.1 )
          
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
std::vector<int> etau_analyzer::getJetCand(int eleIndex, int tauIndex, int ele2Index){
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
      double lepton1Phi=elePhi->at(eleIndex);
      double lepton1Eta= eleEta->at(eleIndex);
      double lepton2Phi=0;double lepton2Eta=0;
      if(ele2Index<0){ lepton2Phi= tau_Phi->at(tauIndex); lepton2Eta=tau_Eta->at(tauIndex); }
      else           { lepton2Phi= elePhi->at(ele2Index); lepton2Eta=eleEta->at(ele2Index); }
      
      double dr_jetEle=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton1Phi, lepton1Eta );
      double dr_jetTau=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton2Phi, lepton2Eta);
      if( dr_jetEle>0.5 && dr_jetTau>0.5 )
	drPassed=true;
	  
      if( (kinematic50 || kinematic30 ) && drPassed==true && jetPUFullID->at(iJet)>>1&1==1)
	tmpCand.push_back(iJet);
    }
  return tmpCand;
}
//The noisy jets are defined as: 20 < pt < 50 && abs(eta) > 2.65 && abs(eta) < 3.139. 
bool etau_analyzer::noisyJet2017(){

  bool noisyJet = false;
  for(int iJet=0;iJet<nJet;iJet++) //Loop over jets
    {
      bool kinematic = false;
      if( jetPt->at(iJet) > 20 
	  && jetPt->at(iJet) < 50
	  && abs(jetEta->at(iJet)) > 2.65
          && abs(jetEta->at(iJet)) < 3.139
	  ) kinematic=true;

      if( kinematic )
        noisyJet=true;
    }
  return noisyJet;
}
int etau_analyzer::thirdLeptonVeto(int eleIndex, int tauIndex){
  // std::vector<int> tmpCand;
  // tmpCand.clear();
  // std::vector<int> tmp2Cand;
  // tmp2Cand.clear();
  // int thirdLepIndex = -1;
  // bool thirdLepVeto=true;
  // for(int iMu=0; iMu < nMu;iMu++)
  //   {
  //     bool kinematic = false;
  //     if( (*muPt)[iMu] > 10.0  && fabs((*muEta)[iMu])< 2.4 && (*muD0)[iMu] < 0.045 && (*muDz)[iMu] < 0.2 ) kinematic = true;
  //     bool muonId =false;
  //     if( muIDbit->at(iMu)>>1&1==1) muonId =true;
  //     bool relative_iso = false;
  //     float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
  //     if( relMuIso < 0.3 ) relative_iso = true;
  //     //if( muIDbit->at(iMu)>>6&1==1) relative_iso =true;
  //     if(muonId==true && kinematic==true && relative_iso==true){
  // 	tmpCand.push_back(iMu);
  //     }                   
  //   }          
  // for(int iEle=0; iEle < nEle;iEle++)
  //   {
  //     bool kinematic = false;
  //     if( (*elePt)[iEle] > 10.0  
  // 	  && fabs((*eleEta)[iEle])< 2.5 
  // 	  && (*eleD0)[iEle] < 0.045 
  // 	  && (*eleDz)[iEle] < 0.2 
  // 	  && eleMissHits->at(iEle) <= 1 && eleConvVeto->at(iEle)==1
  // 	  ) kinematic = true;
  //     bool electronId =false;
  //     if( eleIDbit->at(iEle)>>8&1==1) electronId =true;
  //     bool relative_iso = false;
  //     float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
  //     if( relEleIso < 0.3 ) relative_iso = true;
  //     if(electronId==true && kinematic==true && relative_iso==true){
  // 	tmp2Cand.push_back(iEle);
  //     }                          
  //   }
  // double deltaRe1=0; double deltaRe2=0;
  // double deltaRm1=0; double deltaRm2=0;
  // bool found_3rdele=false; bool found_3rdmu=false;
  // if(tmpCand.size() > 0 )
  //   { 
  //     deltaRm1 = delta_R(elePhi->at(eleIndex),eleEta->at(eleIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
  //     deltaRm2 = delta_R(tau_Phi->at(tauIndex),tau_Eta->at(tauIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
  //     if(deltaRm1>0.5 && deltaRm2>0.5 ){
  // 	thirdLepIndex = tmpCand[0]; found_3rdmu=true;
  //     }
  //   }
  // if(tmp2Cand.size() > 0 )
  //   {
  //     deltaRe1 = delta_R(elePhi->at(eleIndex),eleEta->at(eleIndex), elePhi->at(tmp2Cand[0]),  eleEta->at(tmp2Cand[0]));
  //     deltaRe2 = delta_R(tau_Phi->at(tauIndex),tau_Eta->at(tauIndex), elePhi->at(tmp2Cand[0]),  eleEta->at(tmp2Cand[0]));
  //     if(deltaRe1>0.5 && deltaRe2>0.5 ){
  // 	thirdLepIndex = tmp2Cand[0];found_3rdele=true;
  //     }
  //   }
  // if( found_3rdmu)
  //   return tmpCand[0];
  // else
  //   return -1;

  std::vector<int> tmpCand;
  tmpCand.clear();
  int thirdLepIndex = -1;
  bool thirdLepVeto=true;
  for(int iMu=0; iMu < nMu;iMu++)
    {
      bool kinematic = false;
      if( (*muPt)[iMu] > 10.0  && fabs((*muEta)[iMu])< 2.4 && (*muD0)[iMu] < 0.045 && (*muDz)[iMu] < 0.2 ) kinematic = true;
      bool muonId =false;
      if( muIDbit->at(iMu)>>1&1==1) muonId =true;
      bool relative_iso = false;
      float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
      if( relMuIso < 0.3 ) relative_iso = true;
      if(muonId==true && kinematic==true && relative_iso==true){
	tmpCand.push_back(iMu);
      }                   
    }
  double deltaRm1=0; double deltaRm2=0; bool found_3rdmu=false;
  if(tmpCand.size() > 0 )
    { 
      deltaRm1 = delta_R(elePhi->at(eleIndex),eleEta->at(eleIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      deltaRm2 = delta_R(tau_Phi->at(tauIndex),tau_Eta->at(tauIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      if(deltaRm1>0.5 && deltaRm2>0.5 ){
	found_3rdmu=true;
      }
    }
  
  if(tmpCand.size() > 0 && found_3rdmu==true)
    { thirdLepIndex = tmpCand[0]; thirdLepVeto=false;}
  return thirdLepIndex;

}




bool etau_analyzer::MatchTriggerFilter(int eleIndex, int tauIndex)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool passFilter = true;
  bool eleTriggerFilterMatch=false;
  int nEleTriggerFilterMatch=0;
  bool tauTriggerFilterMatch=false;
  int nTauTriggerFilterMatch=0;
  int nEleDoubleTriggerFilterMatch=0; bool eleDoubleTrgsMatch=false;
  for(int ifilter=39;ifilter<56;ifilter++)
    {
      if(eleFiredSingleTrgs->at(eleIndex)>>ifilter&1==1)
	{
	  eleTriggerFilterMatch=true;
	  nEleTriggerFilterMatch++;
	}
    }
  for(int ifilter=12;ifilter<20;ifilter++)
    {
      if(eleFiredDoubleTrgs->at(eleIndex)>>ifilter&1==1)
        {
          eleDoubleTrgsMatch=true;
          nEleDoubleTriggerFilterMatch++;
        }
    }
  if( (HLTEleMuX>>3&1 == 1 )
      || (HLTEleMuX>>61&1 == 1 )
      || (HLTEleMuX>>5&1 == 1  ) 
      || (HLTTau>>1&1 ==1  )
      ) 
    return true;
  else
    return false;
}



