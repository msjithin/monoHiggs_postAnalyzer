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
  std::vector<int> tauCand;       tauCand.clear();
  
  std::vector<int> higgsCand;     higgsCand.clear();
  std::vector<int> reco_ele;      reco_ele.clear(); 
  std::vector<int> reco_tau;      reco_tau.clear(); 
  
  std::vector<int> eleGenCand;    eleGenCand.clear();
  std::vector<int> muGenCand;     muGenCand.clear();   
  std::vector<int> tauGenCand;    tauGenCand.clear();
  std::vector<int> tauhGenCand;   tauhGenCand.clear();
  std::vector<int> tauNeuGenCand; tauNeuGenCand.clear();
  
  std::vector<int> skimmedElectron;   skimmedElectron.clear(); 
  std::vector<int> skimmedTau;    skimmedTau.clear();
  
  TString sample = TString(SampleName);
  
  int nHiggs = 0;
  int nHToMuTau = 0;
  int found_mt = 0;
  int eleCand_1=0; int eleCand_2=0;int eleCand_3=0;
  int tauCand_1=0; int tauCand_2=0;int tauCand_3=0;
  
  bool fill_hist = false;
  bool isMC = false;
  if( _isMC_=="MC" ) { isMC=true; fill_hist=true; }
  else if ( _isMC_=="DATA" ) { isMC=false; fill_hist=false; }
  
  Double_t  Pt_Bins[26]={0.0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
  Double_t  Pt_Bins_highPt[21]={100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
  
  TH1F* h_cutflow=new TH1F("cutflow", "cutflow", 10, 0, 10); h_cutflow->Sumw2();
  TH1F* h_cutflow_n=new TH1F("cutflow_n", "cutflow_n", 8, 0, 8);h_cutflow_n->Sumw2();
  TH1F* h_cutflow_n_fr=new TH1F("cutflow_n_fr", "cutflow_n_fr", 8, 0, 8);h_cutflow_n_fr->Sumw2();
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
       higgsCand.clear();
       eleCand.clear();
       tauCand.clear();
       reco_ele.clear(); 
       reco_tau.clear();  
       eleGenCand.clear();   
       muGenCand.clear();
       tauGenCand.clear();
       tauhGenCand.clear();
       tauNeuGenCand.clear();
       skimmedElectron.clear();
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

       int report_i=0;
       numberOfEvents+=event_weight;
       event_weight=inspected_event_weight;
       //cout<<jentry<< "  nMC : "<<nMC<<endl;
       //cout<<jentry<< "  genMatching : "<<genMatch2->size()<<endl;
      
       if(isMC)
	 pileup_sf = h_pileup->GetBinContent(h_pileup->GetXaxis()->FindBin(puTrue->at(0)));
       event_weight = event_weight*pileup_sf;
       if( isGoodVtx==false ) continue;
       //if( noisyJet2017()==true ) continue;
       
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
       event_weight=event_weight*kfactor;
       //cout<<" kfactor = "<<kfactor<<endl;
       
       if (isMC) genMatching = gen_matching();
       else genMatching = 0;
       if(debug)cout<<"this worked Line 270"<<endl;
       //// adding skimming
       /////
       int leading_muon = -1; float leading_mPt=0;
       int leading_tau = -1;  float leading_tPt=0;
       ////// reco selection begin
       //if(debug)cout<<"this worked Line 344"<<endl;
       if(metFilters==0)
	 {
	   if(debug)cout<<"metfilters selected"<<endl;
	   if (isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed+=event_weight;
	   if(debug)cout<<"genweight applied"<<endl;
	   if( HLTEleMuX>>5&1==1 
	       //|| HLTEleMuX>>59&1==1
	       || HLTTau>>1&1==1
	       )
	     {
	       nSingleTrgPassed+=event_weight;
	       if(debug)cout<<"trigger selected"<<endl;
	       eleCand = getEleCand(24,2.1);  ///// ele selected 
	       if( eleCand.size() >0 ) 
		 { 
		   nGoodMuonPassed+=event_weight;
		   if(debug)cout<<"this worked Line 296"<<endl;
		   tauCand = getTauCand(30,2.3);
		   if( tauCand.size()>0 ) 
		     {
		       if( HLTTau>>1&1 == 1 )
			 if(! ( tau_Pt->at(tauCand[0]) >32 && abs(tau_Eta->at(tauCand[0])) < 2.1 ))
			   continue;
		       nGoodTauPassed+=event_weight;
		       reco_ele.clear();reco_tau.clear();
		       reco_ele=eleCand; reco_tau=tauCand;
		       if (  eleCharge->at(reco_ele[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
			 {
			   nGoodMuTauPassed+=event_weight;
			   if(debug)cout<<"this worked Line 336"<<endl;
			   afterSF1+=event_weight;
			   
			   if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , false , isMC , debug ) <<endl;
			   if (isMC) event_weight = event_weight * getScaleFactors( reco_ele[0] , reco_tau[0] , false , isMC , debug );
			   
			   afterSF4+=event_weight;
			   if( thirdLeptonVeto() < 0 )
			     {
			       nPassedThirdLepVeto+=event_weight;
			       
			       if( passBjetVeto() == true)
				 {
				   nPassedBjetVeto+=event_weight;
				   
				   double deltaR = delta_R(elePhi->at(reco_ele[0]),eleEta->at(reco_ele[0]), tau_Phi->at(reco_tau[0]),  tau_Eta->at(reco_tau[0]));
				   if(deltaR > 0.5 )
				     {
				       nDeltaRPassed+=event_weight;
				       if(isMC==false)event_weight=1.0;
				       if(debug)cout<<"this worked Line 442"<<endl;
				       fillHist("5", reco_ele[0], reco_tau[0], event_weight);
				       double mT_eleMet = TMass_F((elePt->at(reco_ele[0])),(elePhi->at(reco_ele[0])),pfMET,pfMETPhi  );
				       if( mT_eleMet < 50)
					 {
					   fillHist("6", reco_ele[0], reco_tau[0], event_weight);
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
       ///// fake background region - antiisolated begin
       //event_weight=inspected_event_weight;
       eleCand.clear(); tauCand.clear();
       if(metFilters==0)
	 {
	   if (isMC)fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed_fr+=event_weight;
	   if(  HLTEleMuX>>5&1==1
		//|| HLTEleMuX>>59&1==1
		|| HLTTau>>1&1==1
		)
	     {
	       nSingleTrgPassed_fr+=event_weight;
	       if(debug)cout<<"trigger selected line 368"<<endl;
	       eleCand = getEleCand(24,2.1);  ///// ele selected 
	       if( eleCand.size() >0 ) 
		 { 
		   nGoodMuonPassed_fr+=event_weight;
		   if(debug)cout<<"this worked Line 367"<<endl;
		   tauCand = getAISRTauCand(30,2.3);
		   if( tauCand.size()>0 ) 
		     {
		       if( HLTTau>>1&1 == 1 )
			 if(! ( tau_Pt->at(tauCand[0]) >32 && abs(tau_Eta->at(tauCand[0])) < 2.1 ))
			   continue;
		       nGoodTauPassed_fr+=event_weight;
		       reco_ele.clear();reco_tau.clear();
		       reco_ele=eleCand; reco_tau=tauCand;
		       if (  eleCharge->at(reco_ele[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
			 {
			   nGoodMuTauPassed_fr+=event_weight;
			   
			   if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , true , isMC , debug ) <<endl;
			   if (isMC) event_weight = event_weight * getScaleFactors( reco_ele[0] , reco_tau[0] , true , isMC , debug );
			   
			   if (debug==true ) std::cout<<"event_weight =  "<< event_weight<<" event number = "<<jentry <<std::endl;
			   event_weight = event_weight* getFR(reco_tau[0]);
			   if( thirdLeptonVeto() < 0 )
			     {
			       nPassedThirdLepVeto_fr+=event_weight;
			       
			       if( passBjetVeto() == true)
				 {
				   nPassedBjetVeto_fr+=event_weight;
				   
				   double deltaR = delta_R(elePhi->at(reco_ele[0]),eleEta->at(reco_ele[0]), tau_Phi->at(reco_tau[0]),  tau_Eta->at(reco_tau[0]));
				   if(deltaR > 0.5 )
				     {
				       nDeltaRPassed_fr+=event_weight;
				       if(debug)cout<<"this worked Line 442"<<endl;
				       fillHist("5_fr", reco_ele[0], reco_tau[0], event_weight);
				       double mT_eleMet = TMass_F((elePt->at(reco_ele[0])),(elePhi->at(reco_ele[0])),pfMET,pfMETPhi  );
                                       if( mT_eleMet < 50)
                                         {
                                           fillHist("6_fr", reco_ele[0], reco_tau[0], event_weight);
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
   
   fileName->cd();
   map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
   map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
   for (; iMap1 != jMap1; ++iMap1)
     nplot1(iMap1->first)->Write();
}

void etau_analyzer::BookHistos(const char* file1, const char* file2)
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


void etau_analyzer::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{
  
  //  h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);

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
      if( relEleIso < 0.1 ) relative_iso = true;
      bool trigger = false;
      if( ( HLTEleMuX>>5&1 == 1 && elePt->at(iEle) > 36 &&  eleFiredSingleTrgs->at(iEle)>>13&1==1 )
	  //|| ( HLTEleMuX>>59&1 == 1 && elePt->at(iEle) > 33 )
	  || ( HLTTau>>1&1 == 1 && elePt->at(iEle) > 25  && elePt->at(iEle) < 33 && fabs(eleEta->at(iEle))< 2.1 )
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
      if( tau_Pt->at(iTau) > tauPtCut 
	  && fabs( tau_Eta->at(iTau))< tauEtaCut 
	  && tau_LeadChargedHadron_dz->at(iTau) < 0.2
	  )kinematic = true;
      //if( tau_IDbits->at(iTau)>>16&1==1 ) tauIsolation=true;
      if( tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1 ) tauIsolation=true; 
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byTightDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      //if( tau_IDbits->at(iTau)>>2&1==1 && tau_IDbits->at(iTau)>>7&1==1 )tau_reject=true;
      if( kinematic==true    
	  && decayModeCut==true   
	  && tauIsolation==true 
	  && tau_reject==true   
	  //&& newDecayModeFinding==true
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
      if( tau_Pt->at(iTau) > tauPtCut 
	  && fabs( tau_Eta->at(iTau))< tauEtaCut 
	  && tau_LeadChargedHadron_dz->at(iTau) < 0.2
  	  )kinematic = true;
      if(  tau_byVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 && tau_byMediumDeepTau2017v2p1VSjet->at(iTau)!=1 ) tauIsolation=true;
      //if( tau_IDbits->at(iTau)>>13&1==1 && !(tau_IDbits->at(iTau)>>16&1==1) ) tauIsolation=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byTightDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      if( kinematic==true    
	  && decayModeCut==true   
	  && tauIsolation==true 
	  && tau_reject==true   
	  //&& newDecayModeFinding==true
     	  )
	{
	  tmpCand.push_back(iTau);
    	}                                                           
    }                                                                                       
  return tmpCand;  
}
std::vector<int> etau_analyzer::getJetCand(){
  std::vector<int> tmpCand;  tmpCand.clear();
  for(int iJet=0;iJet<nJet;iJet++) //Loop over jets
    {
      bool kinematic = false;
      bool jet_id = false;
      bool looseID = false;
      if( jetPt->at(iJet) > 30 && abs(jetEta->at(iJet))<4.7) kinematic=true;
      if( jetPt->at(iJet) >= 50 && jetID->at(iJet)>>0&1==1 ) jet_id = true;
      if( jetPt->at(iJet) < 50  && jetPUFullID->at(iJet)>>1&1==1 ) looseID=true;
      
      if( kinematic && (jet_id || looseID) )
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
int etau_analyzer::thirdLeptonVeto(){
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
      //if( muIDbit->at(iMu)>>6&1==1) relative_iso =true;
      if(muonId==true && kinematic==true && relative_iso==true){
	tmpCand.push_back(iMu);
      }                   
    }          
  if(tmpCand.size() > 0){ thirdLepIndex = tmpCand[0]; thirdLepVeto=false;}
  return thirdLepIndex;
  
}
                                                                                    

double etau_analyzer::dR(int mu_index, int tau_index)
{
  double deltaeta = abs(muEta->at(mu_index) - tau_Eta->at(tau_index));
  double muonPhi = muPhi->at(mu_index);
  double tauPhi = tau_Phi->at(tau_index);

  double deltaphi = DeltaPhi(muonPhi, tauPhi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}

double etau_analyzer::delta_R(float phi1, float eta1, float phi2, float eta2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar   = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}



double etau_analyzer::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
}

float etau_analyzer::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
    return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

float etau_analyzer::TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float totalTMass = (a + b+ met).M();
  return totalTMass;
}


float etau_analyzer::VisMass_F(TLorentzVector a, TLorentzVector b){
  float visibleMass = (a + b).M();
  return visibleMass;
}

float etau_analyzer::pTvecsum_F(float pt1, float pt2, float phi1, float phi2) {
  float pt_vecSum = sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2), 2) + pow(pt1*sin(phi1) + pt2*sin(phi2), 2));
  return pt_vecSum;
}

bool etau_analyzer::passBjetVeto()
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
    if( jetPt->at(iJets) > 25  && abs(jetEta->at(iJets)) < 2.4 && jetID->at(iJets)>>1&1==1 ){
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
std::vector<int> etau_analyzer::found_higgs(){ 
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
std::vector<int> etau_analyzer::found_muon(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int i=0; i<nMC;i++){
    if (fabs(mcPID->at(i))==13  )tmpCand.push_back(i); ///&& mcStatus->at(i)==1 
  }
  return tmpCand;
}
std::vector<int> etau_analyzer::found_electron(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int i=0; i<nMC;i++){
    if (fabs(mcPID->at(i))==11 )tmpCand.push_back(i); ///&& mcStatus->at(i)==1
  }
  return tmpCand;
}

std::vector<int> etau_analyzer::found_tau(){
  std::vector<int> tmpCand;    tmpCand.clear(); 
  for(int i=0; i<nMC;i++){
    if ( fabs(mcPID->at(i)) ==15 )tmpCand.push_back(i);
  }
  return tmpCand;
}
std::vector<int> etau_analyzer::found_tauh(){
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

std::vector<int> etau_analyzer::found_tauNeu(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int i=0; i<nMC;i++){
    if (fabs(mcPID->at(i))==16)tmpCand.push_back(i);
  }
  return tmpCand;
}
bool etau_analyzer::skimming_Htt(){
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

int etau_analyzer::gen_matching(){
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
std::vector<int> etau_analyzer::getGenMu(){
  std::vector<int> tmpCand;
  tmpCand.clear();
  int count1=0; int count2=0;
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>2&1==1 || genMatch2->at(imc)>>4&1==1 ) {tmpCand.push_back(imc); count1++;}
    if( (genMatch2->at(imc)>>2&1==1 || genMatch2->at(imc)>>4&1==1) && fabs(mcPID->at(imc))==13 ){count2++;}
  }
  //cout<<"count1:"<<count1<<"  count2:"<<count2<<endl;
  return tmpCand; 
}
float etau_analyzer::exponential(float x,float a,float b,float c) {
  return a * TMath::Exp(-b * x) + c;
}
double etau_analyzer::getScaleFactors( int eleIndex , int tauIndex, bool fakeBkg , bool isMC, bool debug)
{
  double rv_sf=1.0;
  double eleRecoSF_corr=1.0;
  double eleEffSF_corr=1.0;
  double eletrgsf=1.0;
  double sf_Zvtx=1.0;
  double sf_tauidSF_m = 1.0;
  double sf_tauidSF_vvvl = 1.0;
  double sf_tauesSF = 1.0;
  double sf_fakeEle = 1.0; double sf_fakeMu = 1.0;
  double sf_taufesSF = 1.0;
  int genMatchTau = 0;
  if(isMC) genMatchTau = gen_matching();
  
  eleRecoSF_corr=h_eleRecoSF_highpt->GetBinContent(h_eleRecoSF_highpt->GetXaxis()->FindBin(eleSCEta->at(eleIndex)),h_eleRecoSF_highpt->GetYaxis()->FindBin(elePt->at(eleIndex)));
  if (debug==true ) std::cout<<"eleRecoSF_corr =  "<< eleRecoSF_corr<<std::endl;
  eleEffSF_corr=h_eleIDSF->GetBinContent(h_eleIDSF->GetXaxis()->FindBin(eleSCEta->at(eleIndex)),h_eleIDSF->GetYaxis()->FindBin(elePt->at(eleIndex)));
  if (debug==true ) std::cout<<"eleEffSF_corr =  "<< eleEffSF_corr<<std::endl;
  if (debug==true ) std::cout<<"This works line 269 "<<std::endl;
  eletrgsf = EletriggerSF(elePt->at(eleIndex), eleEta->at(eleIndex));
  if (debug==true ) std::cout<<"eletrgsf =  "<< eletrgsf << "Line 271"<<std::endl;
  
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
  if(fakeBkg)
    rv_sf = eleRecoSF_corr * eleEffSF_corr * eletrgsf * sf_Zvtx * sf_tauidSF_m * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF ;
  else
    rv_sf = eleRecoSF_corr * eleEffSF_corr * eletrgsf * sf_Zvtx * sf_tauidSF_m * sf_tauidSF_vvvl * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF ;
  
  return rv_sf;

}
double  etau_analyzer::getFR(int tauIndex){
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

void etau_analyzer::fillHist( string histNumber , int eleIndex, int tauIndex, float event_weight){
  string hNumber = histNumber;
  plotFill("elePt_"+hNumber,  elePt->at(eleIndex) , 40 , 20 , 100,  event_weight);
  plotFill("eleEta_"+hNumber, eleEta->at(eleIndex), 25, -2.5, 2.5,  event_weight);
  plotFill("elePhi_"+hNumber, elePhi->at(eleIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("eleDz_"+hNumber,  eleDz->at(eleIndex), 40, -0.5, 0.5,  event_weight);
  plotFill("eleD0_"+hNumber,  eleD0->at(eleIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("electronID_"+hNumber, eleIDbit->at(eleIndex)>>8&1, 4, -2, 2,  event_weight); // electronID
  float relEleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (elePt->at(eleIndex));
  plotFill("relEleIso_"+hNumber, relEleIso, 15, 0, 0.3,  event_weight);
  plotFill("eleCharge_"+hNumber, eleCharge->at(eleIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 25, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), tau_Phi->at(tauIndex), tau_Eta->at(tauIndex));
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand();
  plotFill("nJet_"+hNumber,  jetCand.size() , 6, 0, 6,  event_weight);
  double HiggsPt = pTvecsum_F(elePt->at(eleIndex),tau_Pt->at(tauIndex), elePhi->at(eleIndex),tau_Phi->at(tauIndex) );
  plotFill("higgsPt_"+hNumber,HiggsPt , 40, 0, 400,  event_weight);
  plotFill("met_"+hNumber, pfMET , 40, 0, 200,  event_weight);
  
  double mT_eleMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
  plotFill("mT_muMet_"+hNumber, mT_eleMet , 30, 0, 150,  event_weight);

  TLorentzVector myTau; 
  myTau.SetPtEtaPhiE(tau_Pt->at(tauIndex),tau_Eta->at(tauIndex),tau_Phi->at(tauIndex), tau_Energy->at(tauIndex));
  TLorentzVector myEle; 
  myEle.SetPtEtaPhiE(elePt->at(eleIndex), eleEta->at(eleIndex), elePhi->at(eleIndex), eleE->at(eleIndex));
  double visMass_mutau = VisMass_F(myTau, myEle);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  double tot_tr_mass = TotTMass_F(myEle, myTau, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}
void etau_analyzer::fillHist( string histNumber , TLorentzVector eleP4, TLorentzVector tauP4, int eleIndex, int tauIndex, float event_weight){
  string hNumber = histNumber;
  plotFill("elePt_"+hNumber,  eleP4.Pt() , 40 , 20 , 100,  event_weight);
  plotFill("eleEta_"+hNumber, eleP4.Eta(), 25, -2.5, 2.5,  event_weight);
  plotFill("elePhi_"+hNumber, eleP4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("eleDz_"+hNumber,  eleDz->at(eleIndex), 40, -0.5, 0.5,  event_weight);
  plotFill("eleD0_"+hNumber,  eleD0->at(eleIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("electronID_"+hNumber, eleIDbit->at(eleIndex)>>8&1, 4, -2, 2,  event_weight); 
  float relEleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (elePt->at(eleIndex));
  plotFill("relEleIso_"+hNumber, relEleIso, 15, 0, 0.3,  event_weight);
  plotFill("eleCharge_"+hNumber, eleCharge->at(eleIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 25, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), tau_Phi->at(tauIndex), tau_Eta->at(tauIndex));
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);

  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand();
  plotFill("nJet_"+hNumber, jetCand.size() , 6, 0, 6,  event_weight);
  double HiggsPt = pTvecsum_F(elePt->at(eleIndex),tau_Pt->at(tauIndex), elePhi->at(eleIndex),tau_Phi->at(tauIndex) );
  plotFill("higgsPt_"+hNumber,HiggsPt , 40, 0, 400,  event_weight);
  plotFill("met_"+hNumber, pfMET , 40, 0, 200,  event_weight);

  double mT_muMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
  plotFill("mT_muMet_"+hNumber, mT_muMet , 30, 0, 150,  event_weight);
  
  double visMass_mutau = VisMass_F(eleP4, tauP4);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  double tot_tr_mass = TotTMass_F(eleP4, tauP4, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
}
float etau_analyzer::EletriggerSF(float pt, float eta){
   float sf = 1.0;
   if(fabs(eta) >= 0.0   && fabs(eta) < 0.8)
     {
       if(pt > 40.0  && pt < 50) sf = 0.79;
       if(pt > 50.0  && pt < 60) sf = 0.82;
       if(pt > 60.0  && pt < 100) sf = 0.85;
       if(pt > 100.0  && pt < 150) sf = 0.87;
       if(pt > 150.0  && pt < 200) sf = 0.88;
       if(pt > 200) sf = 0.89;
     }
   if(fabs(eta) >= 0.8   && fabs(eta) < 1.442 ) 
     {
       if(pt > 40.0  && pt < 50) sf = 0.77;
       if(pt > 50.0  && pt < 60) sf = 0.81;
       if(pt > 60.0  && pt < 100) sf = 0.85;
       if(pt > 100.0  && pt < 150) sf = 0.87;
       if(pt > 150.0  && pt < 300) sf = 0.89;
       if(pt > 300.0) sf = 0.87;
     }
   if(fabs(eta) >= 1.442   && fabs(eta) < 1.56 ) 
     {
       if(pt > 40.0  && pt < 50) sf = 0.73;
       if(pt > 50.0  && pt < 60) sf = 0.75;
       if(pt > 60.0  && pt < 100) sf = 0.76;
       if(pt > 100.0  && pt < 150) sf = 0.72;
       if(pt > 150.0  && pt < 300) sf = 0.78;
       if(pt > 300.0) sf = 0.67;
     }
   if(fabs(eta) >= 1.56   && fabs(eta) < 2.0 ) 
     {
       if(pt > 40.0  && pt < 50) sf = 0.80;
       if(pt > 50.0  && pt < 60) sf = 0.84;
       if(pt > 60.0  && pt < 100) sf = 0.87;
       if(pt > 100.0  && pt < 150) sf = 0.88;
       if(pt > 150.0  && pt < 300) sf = 0.89;
       if(pt > 300.0) sf = 0.87;
     }
   if(fabs(eta) >= 2.0   && fabs(eta) < 2.5 ) 
     {
       if(pt > 40.0  && pt < 50) sf = 0.73;
       if(pt > 50.0  && pt < 60) sf = 0.78;
       if(pt > 60.0  && pt < 100) sf = 0.83;
       if(pt > 100.0  && pt < 150) sf = 0.86;
       if(pt > 150.0  && pt < 300) sf = 0.89;
       if(pt > 300.0) sf = 0.86;
     }
   return sf;
 }
