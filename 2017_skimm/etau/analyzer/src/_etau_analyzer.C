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
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"

#include "commonFunctions.h"
//#include "RecoilCorrector.cc"
#include "ApplyFF.h"
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
  t.Loop(maxEvents,reportEvery, SampleName );
  //delete myMap1;
  cout<<" Outpt written to "<<outputfile<<endl;
  sw.Stop();
  sw.Print();
  return 0;
}

void etau_analyzer::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
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
  int hltele61=0; int eleF45=0;  int eleF53=0;  int eleF54=0;  int eleF55=0;
  int tauF11=0;  int tauF12=0;  int tauF16=0;
  if (fChain == 0) return;
  int genMatching=0; 
  int thirdLeptonIndex=-1;
  std::vector<int> muonGen;
  std::vector<int> eleCand;        eleCand.clear();
  std::vector<int> ele2Cand;       ele2Cand.clear();
  std::vector<int> tauCand;        tauCand.clear();
  std::vector<int> jetCand;       jetCand.clear();
  
  TString sample = TString(SampleName);
  int nHiggs = 0;
  bool fill_hist = false;
  
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
   
   TFile frawff("sf_files/ComputeFF2018/ff_files_et_2017/uncorrected_fakefactors_et.root");
   TF1* ff_qcd_0jet=(TF1*) frawff.Get("rawFF_et_qcd_0jet");
   TF1* ff_qcd_1jet=(TF1*) frawff.Get("rawFF_et_qcd_1jet");
   TF1* ff_w_0jet=(TF1*) frawff.Get("rawFF_et_w_0jet");
   TF1* ff_w_1jet=(TF1*) frawff.Get("rawFF_et_w_1jet");
   TF1* ff_tt_0jet=(TF1*) frawff.Get("mc_rawFF_et_tt");

   TFile fmvisclosure ("sf_files/ComputeFF2018/ff_files_et_2017/FF_corrections_1.root");
   TF1* mvisclosure_qcd=(TF1*) fmvisclosure.Get("closure_mvis_et_qcd");
   TF1* mvisclosure_w=(TF1*) fmvisclosure.Get("closure_mvis_et_w");
   TF1* mvisclosure_tt=(TF1*) fmvisclosure.Get("closure_mvis_et_ttmc");

   TFile fosssclosure ("sf_files/ComputeFF2018/ff_files_et_2017/FF_QCDcorrectionOSSS.root");
   TF1* osssclosure_qcd=(TF1*) fosssclosure.Get("closure_OSSS_mvis_et_qcd");
   TF1* mtclosure_w=(TF1*) fosssclosure.Get("closure_mt_et_w");





   
   Long64_t nentries = fChain->GetEntries();
   if ( is_MC==true ) std::cout<<".... MC file ..... "<<std::endl;
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
       	      
       eleCand.clear();
       tauCand.clear();
       jetCand.clear();
       int reco_ele=0; int reco_ele2=0; int reco_tau=0;
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       double inspected_event_weight = 1.0; 
       if(is_MC)	 fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight/fabs(genWeight) : inspected_event_weight = 0.0;
       nInspected_genWeighted += inspected_event_weight;  
       nInspected += 1; 
       double event_weight=1.0;
       double weight=1.0;
       double applySf=1.0;
       double kfactor = 1.0;
       double nlo_ewk = 1.0;
       double nlo_qcd_binned = 1.0;
       double nlo_qcd = 1.0;
       double nnlo_qcd =1.0;
       int bosonPID;
       double bosonPt=0.0;
       bool Wfound = false;
       double pileup_sf = 1.0;
       bool passSingleTriggerPaths=false;
       bool passCrossTrigger=false;
       int report_i=0;
       bool eleTriggerFilterMatch=false;
       bool tauTriggerFilterMatch=false;
       bool eleTau_selector=false;
       bool eleEle_selector=false;
       bool Ztt_selector=false;

       numberOfEvents+=weight;
       if(is_MC) weight=inspected_event_weight;
       else weight=1.0;
       if(is_MC)
       	 pileup_sf = h_pileup->GetBinContent(h_pileup->GetXaxis()->FindBin(puTrue->at(0)));
       weight = weight*pileup_sf;
       if(is_MC)
       	 weight=weight*prefiringweight;
       if( isGoodVtx==false ) continue;
       

       eleTau_selector=true;
       eleEle_selector=true;
       ////Trigger bit selection
       
       if( ( (HLTEleMuX>>3&1 == 1 )   //HLT_Ele27_WPTight_Gsf_v (HLTEleMuX>>3&1 == 1 )
	     || (HLTEleMuX>>61&1 == 1)  //HLT_Ele32_WPTight_Gsf_v
	     || (HLTEleMuX>>5&1 == 1)   //HLT_Ele35_WPTight_Gsf_v
	     ))
       	 passSingleTriggerPaths=true;  //
       
       if( ( HLTTau>>1&1 == 1 ) )      //HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1
       	 passCrossTrigger=true;
       
       if(is_MC && (found_Wjet_sample || found_DYjet_sample)){
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
       /////
       if(debug)cout<<"entry # : "<<jentry<<endl;
       event_weight=weight;

       if(debug)cout<<"reco selections begin"<<endl;
       eleCand.clear(); ele2Cand.clear();  tauCand.clear();
       ////// reco selection begin
       if(debug)cout<<"signal region DY->ll -  isolated begin"<<endl;
       ////// DY Z-> ll signal region -  isolated begin
       bool dy_ll_genmatching=false;
       //dy_ll_genmatching=true;
       TLorentzVector dau2;  TLorentzVector dau1;
       TLorentzVector initialMet; TLorentzVector initialMetTau; TLorentzVector corrMet;
       if(metFilters==0 )
	 {
	   if(debug)cout<<"metfilters selected"<<endl;
	   if(is_MC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed_dyll+=event_weight;
	   makeTestPlot("a_dyll", 0,0,0,event_weight);
	   if(debug)cout<<"genweight applied"<<endl;
	   if(   passSingleTriggerPaths || passCrossTrigger )
	     {
	       nSingleTrgPassed_dyll+=event_weight;
	       if(debug)cout<<"trigger selected"<<endl;
	       makeTestPlot("b_dyll", 0,0,0,event_weight);
	       eleCand = getEleCand(25.0,2.1);  ///// ele selected
	       if( eleCand.size() >0 ) 
		 { 
		   nGoodMuonPassed_dyll+=event_weight;
		   if(debug)cout<<"this worked Line 443"<<endl;
		   makeTestPlot("c_dyll", 0,0,0,event_weight);
		   
		   tauCand = getTauCand(30,2.3);
		   if( tauCand.size()>0  ) 
		     {
		       nGoodTauPassed_dyll+=event_weight;
		       if(debug)cout<<"this worked Line 424"<<endl;
		       makeTestPlot("d_dyll", 0,0,0,event_weight);
		       
		       if( passDiElectronVeto(eleCand[0])==true 
			   && eVetoZTTp001dxyz(eleCand[0], tauCand[0])
			   && mVetoZTTp001dxyz(eleCand[0], tauCand[0])
			   ) Ztt_selector=true;
		       else Ztt_selector=false;
		     
		   
		       if(Ztt_selector) 
			 {
			   reco_ele=eleCand[0]; reco_tau=tauCand[0];
			   makeTestPlot("e_dyll", 0,0,0,event_weight);
			   if( is_MC )
			     {
			       if(  myGenMaching(reco_tau)<5 &&  myGenMaching1(reco_ele)<5 )
				 dy_ll_genmatching=true;
			       else
				 dy_ll_genmatching=false;
			     }
			   else if(!is_MC)
			     dy_ll_genmatching=true;
		       
			   if ( eleCharge->at(reco_ele) * tau_Charge->at(reco_tau) < 0  &&  dy_ll_genmatching ) 
			     {
			       nGoodMuTauPassed_dyll+=event_weight;
			   
			   
			       if ( MatchTriggerFilter(reco_ele, reco_tau) )
				 {
			       
				   if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
				   if(is_MC)
				     initialMet=MetRecoilCorrections(reco_ele, reco_tau);
				   else
				     initialMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
			       
				   
				   dau2.SetPtEtaPhiE(tau_Pt->at(reco_tau),tau_Eta->at(reco_tau),tau_Phi->at(reco_tau), tau_Energy->at(reco_tau));
				   dau1.SetPtEtaPhiE(elePt->at(reco_ele), eleEta->at(reco_ele), elePhi->at(reco_ele), eleE->at(reco_ele));
				   initialMetTau=initialMet+dau2;
				   
				   if(is_MC)
				     applyESCorrections(dau1,dau2, reco_ele, reco_tau, false,
							dau1, dau2);
				   
				   corrMet = initialMetTau-dau2;
				   if(!is_MC)
				     corrMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);

				   applySf=1.0;
				   if(is_MC)
				     applySf=  getScaleFactors( elePt->at(reco_ele),
								tau_Pt->at(reco_tau),
								eleEta->at(reco_ele),
								tau_Eta->at(reco_tau),
								tau_DecayMode->at(reco_tau),
								myGenMaching(reco_tau),
								false
								);
			       
				   // if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , false , is_MC , debug ) <<endl;
				   event_weight = event_weight * applySf;
				   afterSF4+=event_weight;
				   makeTestPlot("f_dyll", 0,0,0,event_weight);
				   if( thirdLeptonVeto(reco_ele , reco_tau)  )
				     {
				       nPassedThirdLepVeto_dyll+=event_weight;
				       makeTestPlot("g_dyll", 0,0,0,event_weight);
				       if( passBjetVetoM(reco_ele , reco_tau)
					   && passBjetVetoL(reco_ele , reco_tau))
					 {
					   nPassedBjetVeto_dyll+=event_weight;
					   makeTestPlot("h_dyll", 0,0,0,event_weight);
					   double deltaR = delta_R(elePhi->at(reco_ele),eleEta->at(reco_ele), tau_Phi->at(reco_tau),  tau_Eta->at(reco_tau));
					   if(deltaR > 0.5 )
					     {
					       nDeltaRPassed_dyll+=event_weight;
					       if(is_MC==false)event_weight=1.0;
					       makeTestPlot("i_dyll", 0,0,0,event_weight);
					       if(debug)cout<<"this worked Line 374"<<endl;
					       fillHist("5_dyll",  reco_ele, reco_tau, false, event_weight);
					   
					   
					       double mT_eleMet = TMass_F( dau1.Pt(), dau1.Phi(),corrMet.Pt(), corrMet.Phi()  );
					       if( mT_eleMet < 50 )
						 {
						   fillHist("6_dyll", reco_ele, reco_tau, false, event_weight);
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
     
       
       /////// signal region DY Z->ll end
       ////// fake background region - antiisolated  DY Z->ll begin
       if(debug)cout<<"moving to fake bkg DY->ll "<<endl;
       event_weight=weight;
       eleCand.clear(); ele2Cand.clear();  tauCand.clear();
       if(metFilters==0 )
	 {
	   if(is_MC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed_dyll_fr+=event_weight;
	   if(  passSingleTriggerPaths || passCrossTrigger  )
	     {
	       nSingleTrgPassed_dyll_fr+=event_weight;
	       if(debug)cout<<"trigger selected"<<endl;
	       eleCand = getEleCand(25,2.1);  ///// ele selected
	       if( eleCand.size() >0 ) 
		 { 
		   nGoodMuonPassed_dyll_fr+=event_weight;
		   tauCand = getAISRTauCand(30,2.3);
		   if( tauCand.size()>0 ) 
		     {
		       nGoodTauPassed_dyll_fr+=event_weight;
		       if(debug)cout<<"fr tau selection passed"<<endl;
			   
		       if( passDiElectronVeto(eleCand[0])==true 
			   && eVetoZTTp001dxyz(eleCand[0], tauCand[0])
			   && mVetoZTTp001dxyz(eleCand[0], tauCand[0])
			   ) Ztt_selector=true;
		       else Ztt_selector=false;
			   
			   
		       if(Ztt_selector) 
			 {
			   reco_ele=eleCand[0]; reco_tau=tauCand[0];
			   if( is_MC )
			     {
			       if(  myGenMaching(reco_tau)<5 &&  myGenMaching1(reco_ele)<5 )
				 dy_ll_genmatching=true;
			       else
				 dy_ll_genmatching=false;
			     }
			   else if(!is_MC)
			     dy_ll_genmatching=true;
			       
			   if ( eleCharge->at(reco_ele) * tau_Charge->at(reco_tau) < 0   &&  dy_ll_genmatching ) 
			     {
			       nGoodMuTauPassed_dyll_fr+=event_weight;
				   
			       if ( MatchTriggerFilter(reco_ele, reco_tau) )
				 {
				       
				       
				       
				   event_weight = event_weight* getFR(reco_tau);
				   if(is_MC)
				     initialMet=MetRecoilCorrections(reco_ele, reco_tau);
				   else
				     initialMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
				       
				   dau2.SetPtEtaPhiE(tau_Pt->at(reco_tau),tau_Eta->at(reco_tau),tau_Phi->at(reco_tau), tau_Energy->at(reco_tau));
				   dau1.SetPtEtaPhiE(elePt->at(reco_ele), eleEta->at(reco_ele), elePhi->at(reco_ele), eleE->at(reco_ele));
				   initialMetTau=initialMet+dau2;
				       
				       
				   if(is_MC)
                                     applyESCorrections(dau1,dau2, reco_ele, reco_tau, true,
							dau1, dau2);		
				   
				   corrMet = initialMetTau-dau2;
				   if(!is_MC)
                                     corrMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);

				   applySf=1.0;
				   if(is_MC)
				     applySf=  getScaleFactors( elePt->at(reco_ele),
                                                                tau_Pt->at(reco_tau),
								eleEta->at(reco_ele),
								tau_Eta->at(reco_tau),
								tau_DecayMode->at(reco_tau),
								myGenMaching(reco_tau),
								true
								);
				       
				   // if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , false , is_MC , debug ) <<endl;
				   event_weight = event_weight * applySf;
				   /////
				   if( thirdLeptonVeto(reco_ele , reco_tau)  )
				     {
				       nPassedThirdLepVeto_dyll_fr+=event_weight;
					   
				       if( passBjetVetoM(reco_ele , reco_tau)
					   && passBjetVetoL(reco_ele , reco_tau) )
					 {
					   nPassedBjetVeto_dyll_fr+=event_weight;
					       
					   double deltaR = delta_R(elePhi->at(reco_ele),eleEta->at(reco_ele), tau_Phi->at(reco_tau),  tau_Eta->at(reco_tau));
					   if(deltaR > 0.5 )
					     {
					       nDeltaRPassed_dyll_fr+=event_weight;
					       if(debug)cout<<"this worked Line 425"<<endl;
					       if(debug)cout<<"evnt_weight = "<<event_weight<<endl;
					       fillHist("5_dyll_fr", reco_ele, reco_tau, true, event_weight);
					       double mT_eleMet = TMass_F( dau1.Pt(), dau1.Phi(),corrMet.Pt(), corrMet.Phi() );
					       if( mT_eleMet < 50 ) 
						 {
						   fillHist("6_dyll_fr", reco_ele, reco_tau, true, event_weight);
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
       if(debug)cout<<"signal region -  isolated begin L523"<<endl;       
       
       Ztt_selector=false;
       bool dy_ltau_genmatching=false;
       ////// signal region -  isolated begin
       if(is_MC)
	 event_weight=weight;
       else
	 event_weight=1.0;
       eleCand.clear(); ele2Cand.clear();  tauCand.clear();
	   
       if(metFilters==0)
	 {
	   if(debug)cout<<"metfilters selected"<<endl;
	   if (is_MC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed+=event_weight;
	   makeTestPlot("a", 0,0,0,event_weight);
	   if(debug)cout<<"genweight applied"<<endl;
	   if( passSingleTriggerPaths || passCrossTrigger  )
	     {
	       nSingleTrgPassed+=event_weight;
	       if(debug)cout<<"trigger selected"<<endl;
	       makeTestPlot("b", 0,0,0,event_weight);
	       eleCand = getEleCand(25.0,2.1);  ///// ele selected 
	       if( eleCand.size() >0 ) 
		 { 
		   nGoodMuonPassed+=event_weight;
		   if(debug)cout<<"this worked Line 526"<<endl;
		       
		   makeTestPlot("c", 0,0,0,event_weight);
		   tauCand = getTauCand(30.0,2.3);
		   if( tauCand.size() >0 )
		     {
		       nGoodTauPassed+=event_weight;
			   
		       reco_ele=eleCand[0]; reco_tau=tauCand[0];
		       makeTestPlot("d", 0,0,0,event_weight);
			   
		       if( passDiElectronVeto(reco_ele)==true
			   && (eVetoZTTp001dxyz(reco_ele, reco_tau))
			   && (mVetoZTTp001dxyz(reco_ele, reco_tau))
			   ) Ztt_selector=true;
		       else Ztt_selector=false;
			   
			   
			   
		       if(Ztt_selector) 
			 {
			   if( is_MC )
			     {
			       if(  myGenMaching(reco_tau)>4 &&  myGenMaching1(reco_ele)<5 )
				 dy_ltau_genmatching=true;
			       else
				 dy_ltau_genmatching=false;
			     }
			   else if(!is_MC)
			     dy_ltau_genmatching=true;
			       
			       
			   if (  eleCharge->at(reco_ele) * tau_Charge->at(reco_tau) < 0   && dy_ltau_genmatching) 
			     {
			       nGoodMuTauPassed+=event_weight;
			       
			       if(debug)cout<<"this worked Line 538"<<endl;
			       afterSF1+=event_weight;
			       makeTestPlot("e", 0,0,0,event_weight);
			       if ( MatchTriggerFilter(reco_ele, reco_tau) )
				 {
				   if(debug)cout<<"this worked Line 534"<<endl;
				   if(is_MC)
				     initialMet=MetRecoilCorrections(reco_ele, reco_tau);
				   else
				     initialMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
				   
				   dau2.SetPtEtaPhiE(tau_Pt->at(reco_tau),tau_Eta->at(reco_tau),tau_Phi->at(reco_tau), tau_Energy->at(reco_tau));
				   dau1.SetPtEtaPhiE(elePt->at(reco_ele), eleEta->at(reco_ele), elePhi->at(reco_ele), eleE->at(reco_ele));
				   initialMetTau=initialMet+dau2;

                                   if(is_MC)
                                     applyESCorrections(dau1,dau2, reco_ele, reco_tau, false, 
                                                        dau1, dau2);

                                   corrMet = initialMetTau-dau2;
				   if(!is_MC)
                                     corrMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);

				   applySf=1.0;
				   if(is_MC)
				     applySf=  getScaleFactors( elePt->at(reco_ele),
                                                                tau_Pt->at(reco_tau),
								eleEta->at(reco_ele),
								tau_Eta->at(reco_tau),
								tau_DecayMode->at(reco_tau),
								myGenMaching(reco_tau),
								false
								);
				   // if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , false , is_MC , debug ) <<endl;
                                   event_weight = event_weight * applySf;
				   afterSF4+=event_weight;
				   makeTestPlot("f", 0,0,0,event_weight);
				   if( thirdLeptonVeto(reco_ele , reco_tau)  )
				     {
				       nPassedThirdLepVeto+=event_weight;
				       makeTestPlot("g", 0,0,0,event_weight);
				       if( passBjetVetoM(reco_ele , reco_tau)
					   && passBjetVetoL(reco_ele , reco_tau) )
					 {
					   nPassedBjetVeto+=event_weight;
					   makeTestPlot("h", 0,0,0,event_weight);
					   if(tau_DecayMode->at(reco_tau)==5 || tau_DecayMode->at(reco_tau)==6) continue;
					   
					   double deltaR = delta_R(elePhi->at(reco_ele),eleEta->at(reco_ele), tau_Phi->at(reco_tau),  tau_Eta->at(reco_tau));
					   //if(Overlap_2(eleP4, tauP4)==false)
					   if(deltaR > 0.5 )
					     {
					       nDeltaRPassed+=event_weight;
					       if(is_MC==false)event_weight=1.0;
					       if(debug)cout<<"this worked Line 558"<<endl;
					       fillHist("5", reco_ele, reco_tau, false, event_weight);
					       makeTestPlot("i", 0,0,0,event_weight);
					       double mT_eleMet = TMass_F( dau1.Pt(), dau1.Phi() , corrMet.Pt(), corrMet.Phi() );
					       if( mT_eleMet < 50 )
						 {
						   fillHist("6", reco_ele, reco_tau, false, event_weight);
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
       //////// signal region end
       if(debug)cout<<"fake background region - antiisolated begin 625"<<endl;
       ///// fake background region - antiisolated begin
       event_weight=weight;
       eleCand.clear(); tauCand.clear();
       if(metFilters==0)
	 {
	   if (is_MC)fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed_fr+=event_weight;
	   makeTestPlot("a_fr", 0,0,0,event_weight);
	   if(  passSingleTriggerPaths || passCrossTrigger  )
	     {
	       nSingleTrgPassed_fr+=event_weight;
	       if(debug)cout<<"trigger selected line 636"<<endl;
	       makeTestPlot("b_fr", 0,0,0,event_weight);
	       eleCand = getEleCand(25.0,2.1);  ///// ele selected 
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

		       reco_ele=eleCand[0]; reco_tau=tauCand[0];
			   
		       if( passDiElectronVeto(reco_ele)==true
			   && (eVetoZTTp001dxyz(reco_ele, reco_tau))
			   && (mVetoZTTp001dxyz(reco_ele, reco_tau))
			   ) Ztt_selector=true;
		       else Ztt_selector=false;
			   
		       if(Ztt_selector) 
			 {
			   // if( is_MC )
			   //   {
			   //     if(  myGenMaching(reco_tau)>4 &&  myGenMaching1(reco_ele)<5 )
			   // 	 dy_ltau_genmatching=true;
			   //     else
			   // 	 dy_ltau_genmatching=false;
			   //   }
			   // else if(!is_MC)
			   dy_ltau_genmatching=true;
			   if (  eleCharge->at(reco_ele) * tau_Charge->at(reco_tau) < 0  && dy_ltau_genmatching) 
			     {
			       nGoodMuTauPassed_fr+=event_weight;
			       makeTestPlot("e_fr", 0,0,0,event_weight);
			       if ( MatchTriggerFilter(reco_ele, reco_tau) )
				 {
				   if(is_MC)
				     initialMet=MetRecoilCorrections(reco_ele, reco_tau);
				   else
				     initialMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);

				   dau2.SetPtEtaPhiE(tau_Pt->at(reco_tau),tau_Eta->at(reco_tau),tau_Phi->at(reco_tau), tau_Energy->at(reco_tau));
				   dau1.SetPtEtaPhiE(elePt->at(reco_ele), eleEta->at(reco_ele), elePhi->at(reco_ele), eleE->at(reco_ele));
				   
				   initialMetTau=initialMet+dau2;

                                   if(is_MC)
                                     applyESCorrections(dau1,dau2, reco_ele, reco_tau, true,
                                                        dau1, dau2);

                                   corrMet = initialMetTau-dau2;
				   if(!is_MC)
                                     corrMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);

				   applySf=1.0;
				   if(is_MC)
				     applySf=  getScaleFactors( elePt->at(reco_ele),
                                                                tau_Pt->at(reco_tau),
								eleEta->at(reco_ele),
								tau_Eta->at(reco_tau),
								tau_DecayMode->at(reco_tau),
								myGenMaching(reco_tau),
								true
								);
				   
                                   // if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , false , is_MC , debug ) <<endl;
                                   event_weight = event_weight * applySf;
				   if (debug==true ) std::cout<<"event_weight =  "<< event_weight<<" event number = "<<jentry <<std::endl;
				   //event_weight = event_weight* getFR(reco_tau);
				   float mt=TMass_F(dau1.Pt(),dau1.Phi(),corrMet.Pt(), corrMet.Phi());
				   float mvis=VisMass_F(dau1, dau2);
				   std::vector<int> jetCand;       jetCand.clear();
				   jetCand=getJetCand(reco_ele, reco_tau);
				   int njets=jetCand.size();
				   //cout<<"njets="<<njets<<endl;
				   float frac_tt=0.01; float frac_qcd=0.24; float frac_w=0.75; 
				   if(njets==0 && corrMet.Pt()>60.0)
				     { frac_w=1.0; frac_qcd=0.0; frac_tt=0.0; }
				   float my_fakefactor = get_ff( dau2.Pt(), mt, mvis, njets, frac_tt, frac_qcd, frac_w, ff_qcd_0jet, ff_qcd_1jet, ff_w_0jet, ff_w_1jet, ff_tt_0jet, mvisclosure_qcd, mvisclosure_w, mvisclosure_tt, mtclosure_w, osssclosure_qcd);
				   //cout<<" smhtt_fakefactor="<<my_fakefactor<<endl;
				   event_weight = event_weight*my_fakefactor;
				   makeTestPlot("f_fr", 0,0,0,event_weight);
				   if( thirdLeptonVeto(reco_ele , reco_tau) )
				     {
				       nPassedThirdLepVeto_fr+=event_weight;
				       makeTestPlot("g_fr", 0,0,0,event_weight);
				       if( passBjetVetoM(reco_ele , reco_tau)
					   && passBjetVetoL(reco_ele , reco_tau) )
					 {
					   nPassedBjetVeto_fr+=event_weight;
					   makeTestPlot("h_fr", 0,0,0,event_weight);
					   if(tau_DecayMode->at(reco_tau)==5 || tau_DecayMode->at(reco_tau)==6) continue;

					   double deltaR = delta_R(elePhi->at(reco_ele),eleEta->at(reco_ele), tau_Phi->at(reco_tau),  tau_Eta->at(reco_tau));
					   if(deltaR > 0.5 )
					     {
					       nDeltaRPassed_fr+=event_weight;
					       makeTestPlot("i_fr", 0,0,0,event_weight);
					       if(debug)cout<<"this worked Line 442"<<endl;
					       fillHist("5_fr", reco_ele, reco_tau, true, event_weight);
					       double mT_eleMet = TMass_F( dau1.Pt(), dau1.Phi(), corrMet.Pt(), corrMet.Phi()  );
					       if( mT_eleMet < 50 )
						 {
						   fillHist("6_fr", reco_ele, reco_tau, true, event_weight);
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
   
   // std::cout<<"HLTTau>>1&1 ==1 = "<<hltele61<<std::endl;
   // std::cout<<"eleF45 = "<<eleF45<<std::endl;
   // std::cout<<"eleF53 = "<<eleF53<<std::endl;
   // std::cout<<"eleF54 = "<<eleF54<<std::endl;
   // std::cout<<"eleF55 = "<<eleF55<<std::endl;
   // std::cout<<"tauF11 = "<<tauF11<<std::endl;
   // std::cout<<"tauF12 = "<<tauF12<<std::endl;
   // std::cout<<"tauF16 = "<<tauF16<<std::endl;

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
  TLorentzVector dau1;
  //Loop over electrons                                                                     
  for(int iEle=0;iEle<nEle;iEle++)
    {
      dau1.SetPtEtaPhiE(elePt->at(iEle), eleEta->at(iEle), elePhi->at(iEle), eleE->at(iEle));
      dau1 = dau1*(eleCalibE->at(iEle)/dau1.E());
      
      bool kinematic = false;
      if( dau1.Pt() > elePtCut  
	  && fabs(dau1.Eta())< eleEtaCut 
	  && fabs(eleD0->at(iEle)) < 0.045
	  && fabs(eleDz->at(iEle)) < 0.2
	  && eleMissHits->at(iEle) <= 1 && eleConvVeto->at(iEle)==1
	  ) kinematic = true;
      bool electronId =false;
      if( eleIDbit->at(iEle)>>8&1==1) electronId =true;
      bool relative_iso = false;    
      float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / dau1.Pt();
      if( relEleIso < 0.15 ) relative_iso = true;
      bool trigger = false;
      if( ( HLTEleMuX>>3&1 == 1 && dau1.Pt() > 28.0  ) 
	  || ( HLTEleMuX>>61&1 == 1 && dau1.Pt() > 33.0  )
	  || ( HLTEleMuX>>5&1 == 1 && dau1.Pt() > 36.0 )
	  || ( HLTTau>>1&1 == 1 && dau1.Pt() > 25.0  && dau1.Pt() < 28.0 && fabs(dau1.Eta())< 2.1 )
	  
	  ) trigger = true;
      if( kinematic && electronId && relative_iso && trigger ){
	tmpCand.push_back(iEle);
      }	
    }                           
  return tmpCand;
  
}
std::vector<int> etau_analyzer::getEle2Cand(double elePtCut, double eleEtaCut, int ele1Index){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over electrons                                                                     
  for(int iEle=0;iEle<nEle;iEle++)
    {
      if(iEle==ele1Index) continue;
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
  TLorentzVector dau2;
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
	  && fabs(tau_Charge->at(iTau))==1
	  )kinematic = true;
      if( tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1 ) tauIsolation=true; 
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byTightDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      if( ( HLTEleMuX>>3&1 == 1 ) 
	  || ( HLTEleMuX>>61&1 == 1 )
	  || ( HLTEleMuX>>5&1 == 1 )
	  || ( HLTTau>>1&1 == 1 && tau_Pt->at(iTau) >35.0 && fabs(tau_Eta->at(iTau)) < 2.1 )
	 
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
      if(  tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 && !(tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1) ) tauIsolation=true;
      //if( tau_IDbits->at(iTau)>>13&1==1 && !(tau_IDbits->at(iTau)>>16&1==1) ) tauIsolation=true;
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
int etau_analyzer::getZCand()
{
  //std::vector<int> tmpCand;  tmpCand.clear();
  for(int iMC=0;iMC<nMC;iMC++) //Loop over mc
    {
      //if(fabs(mcPID->at(iMC))==23 && mcStatus->at(iMC)==62) tmpCand.push_back(iMC); 
      //if(fabs(mcPID->at(iMC))==23 ) tmpCand.push_back(iMC);
      if(fabs(mcPID->at(iMC))==23 && mcStatus->at(iMC)==62)
	{
	  return iMC;
	}
    }
  return -1;
}
std::vector<int> etau_analyzer::getJetCand(int eleIndex, int tauIndex){
  std::vector<int> tmpCand;  tmpCand.clear();
  for(int iJet=0;iJet<nJet;iJet++) //Loop over jets
    {
      bool kinematic30 = false; bool foundNoisyJets=false;
      bool kinematic50 = false; bool passLoosePUID=false;
      bool jet_id = false; bool drPassed=false;
      if( jetPt->at(iJet) > 30
          && abs(jetEta->at(iJet))<4.7
          && (jetID->at(iJet)>>0&1)==1
          //&& jetPUFullID->at(iJet)>>1&1==1
          ) kinematic30=true;
      if(jetPt->at(iJet) < 50
         && abs(jetEta->at(iJet))>2.65
         && abs(jetEta->at(iJet))<3.139
         && (jetID->at(iJet)>>0&1)==1
         ) foundNoisyJets=true;

      // if( jetPt->at(iJet) > 50
      //          && abs(jetEta->at(iJet))<2.65
      //          && abs(jetEta->at(iJet))>3.139
      //          && (jetID->at(iJet)>>0&1)==1
      //          //&& jetPUFullID->at(iJet)>>1&1==1
      //          ) kinematic50=true;
      // else if( jetPt->at(iJet) > 30
      //               && abs(jetEta->at(iJet))<4.7
      //               && (jetID->at(iJet)>>0&1)==1
      //               //&& jetPUFullID->at(iJet)>>1&1==1
      //               ) kinematic30=true;
      if( jetPt->at(iJet) < 50 )
        {
          if(jetPUFullID->at(iJet)>>1&1==1 )
            passLoosePUID=true;
          else
            passLoosePUID=false;
        }
      else if (jetPt->at(iJet) > 50 )
        passLoosePUID=true;

      double lepton1Phi=elePhi->at(eleIndex);
      double lepton1Eta= eleEta->at(eleIndex);
      double lepton2Phi=0;double lepton2Eta=0;
      lepton2Phi= tau_Phi->at(tauIndex); lepton2Eta=tau_Eta->at(tauIndex);
      double dr_jetEle=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton1Phi, lepton1Eta );
      double dr_jetTau=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton2Phi, lepton2Eta);
      if( dr_jetEle>0.5 && dr_jetTau>0.5 )
        drPassed=true;

      //if( (kinematic50 || kinematic30 ) && drPassed==true  )
      if(kinematic30 && !foundNoisyJets && passLoosePUID && drPassed==true)
        tmpCand.push_back(iJet);
    }
  return tmpCand;
}
//The noisy jets are defined as: 20 < pt < 50 && abs(eta) > 2.65 && abs(eta) < 3.139. 
bool etau_analyzer::noisyJet2017(int reco_ele, int reco_tau){
  
std::vector<int> jetCand;       jetCand.clear();
 jetCand=getJetCand(reco_ele, reco_tau);
 if(jetCand.size()>0)
   {
     int iJet=jetCand[0];
     if( jetRawPt->at(iJet) > 20 
	 && jetRawPt->at(iJet) < 50
	 && abs(jetEta->at(iJet)) > 2.65
	 && abs(jetEta->at(iJet)) < 3.139
	 )
       return true;
   }
 return false;
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
      if(muonId==true && kinematic==true && relative_iso==true){
	tmpCand.push_back(iMu);
      }                   
    }          
  if(tmpCand.size() > 0){ thirdLepIndex = tmpCand[0]; thirdLepVeto=false;}
  return thirdLepIndex;
  
}
bool etau_analyzer::thirdLeptonVeto(int eleIndex, int tauIndex){
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
	return false;
      }
    }
  else
    return true;
  
}
bool etau_analyzer::thirdLeptonVeto(int eleIndex, int tauIndex, int ele2Index){
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
  double deltaRm1=0; double deltaRm2=0; double deltaRm3=0;
  if(tmpCand.size() > 0 )
    { 
      deltaRm1 = delta_R(elePhi->at(eleIndex),eleEta->at(eleIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      deltaRm2 = delta_R(tau_Phi->at(tauIndex),tau_Eta->at(tauIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      deltaRm3 = delta_R(elePhi->at(ele2Index),eleEta->at(ele2Index), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      if(deltaRm1>0.5 && deltaRm2>0.5 && deltaRm3>0.5){
	return false;
      }
    }
  else
    return true;
  
}
                                                                                    

// double etau_analyzer::dR(int ele_index, int tau_index)
// {
//   double deltaeta = abs(eleEta->at(ele_index) - tau_Eta->at(tau_index));
//   double electronPhi = elePhi->at(ele_index);
//   double tauPhi = tau_Phi->at(tau_index);

//   double deltaphi = DeltaPhi(electronPhi, tauPhi);
//   double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
//   return deltar;
  
// }

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
  //return  sqrt(2.0*LepPt*met*(1.0-cos(DeltaPhi(LepPhi, metPhi))));
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
float etau_analyzer::pTvecsum_F(TLorentzVector a, TLorentzVector b, TLorentzVector c) {
  float pt_vecSum = (a + b+ c).Pt();
  return pt_vecSum;
}

bool etau_analyzer::passBjetVeto(int eleIndex, int tauIndex)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  double lepton1Phi=elePhi->at(eleIndex);
  double lepton1Eta= eleEta->at(eleIndex);
  double lepton2Phi= tau_Phi->at(tauIndex); double lepton2Eta=tau_Eta->at(tauIndex);   
  double dr_jetEle=0.0; double dr_jetTau=0.0; 
  
  for(int iJets=0; iJets<nJet ; iJets++){
    bool particles_separated=false;
    dr_jetEle =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton1Phi, lepton1Eta );
    dr_jetTau =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton2Phi, lepton2Eta);
    if( dr_jetEle>0.5 && dr_jetTau>0.5) { particles_separated=true;}
      
    if( jetPt->at(iJets) > 25  
	&& abs(jetEta->at(iJets)) < 2.4 
	&& jetID->at(iJets)>>0&1==1 
	&& jetPUFullID->at(iJets)>>1&1==1
	&& particles_separated==true
	){
      tmpJetCand.push_back(iJets);
    }
  }
  if(tmpJetCand.size() ==1 ){
    // atleast one jet ==> events pass medium 
    if( (jetDeepCSVTags_b->at(tmpJetCand[0]) + jetDeepCSVTags_bb->at(tmpJetCand[0])) > 0.4941  )
      return veto = false;
  }
  else if(tmpJetCand.size() > 1){
    // atleast 2 jets ==> events pass loose
    if( (jetDeepCSVTags_b->at(tmpJetCand[0]) + jetDeepCSVTags_bb->at(tmpJetCand[0])) > 0.1522   )
      return veto = false;
  }
  return veto=true;
}
bool etau_analyzer::passBjetVetoM(int eleIndex, int tauIndex)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  double lepton1Phi=elePhi->at(eleIndex);
  double lepton1Eta= eleEta->at(eleIndex);
  double lepton2Phi= tau_Phi->at(tauIndex); double lepton2Eta=tau_Eta->at(tauIndex);   
  double dr_jetEle=0.0; double dr_jetTau=0.0; 
  
  for(int iJets=0; iJets<nJet ; iJets++){
    bool particles_separated=false;
    dr_jetEle =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton1Phi, lepton1Eta );
    dr_jetTau =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton2Phi, lepton2Eta);
    if( dr_jetEle>0.5 && dr_jetTau>0.5) { particles_separated=true;}
      
    if( jetPt->at(iJets) > 20
	&& abs(jetEta->at(iJets)) < 2.4 
	&& jetID->at(iJets)>>0&1==1 
	&& jetPUFullID->at(iJets)>>1&1==1
	&& particles_separated==true
	&& (jetDeepCSVTags_b->at(iJets) + jetDeepCSVTags_bb->at(iJets)) > 0.4941
	){
      tmpJetCand.push_back(iJets);
    }
  }
  if(tmpJetCand.size() >0 ){
    return veto=false;
  }
  return veto=true;
}
bool etau_analyzer::passBjetVetoL(int eleIndex, int tauIndex)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  double lepton1Phi=elePhi->at(eleIndex);
  double lepton1Eta= eleEta->at(eleIndex);
  double lepton2Phi= tau_Phi->at(tauIndex); double lepton2Eta=tau_Eta->at(tauIndex);   
  double dr_jetEle=0.0; double dr_jetTau=0.0; 
  
  for(int iJets=0; iJets<nJet ; iJets++){
    bool particles_separated=false;
    dr_jetEle =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton1Phi, lepton1Eta );
    dr_jetTau =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton2Phi, lepton2Eta);
    if( dr_jetEle>0.5 && dr_jetTau>0.5) { particles_separated=true;}
      
    if( jetPt->at(iJets) > 20  
	&& abs(jetEta->at(iJets)) < 2.4 
	&& jetID->at(iJets)>>0&1==1 
	&& jetPUFullID->at(iJets)>>1&1==1
	&& particles_separated==true
	&& (jetDeepCSVTags_b->at(iJets) + jetDeepCSVTags_bb->at(iJets)) > 0.1522 
	){
      tmpJetCand.push_back(iJets);
    }
  }
  if(tmpJetCand.size() > 1 ){
    // atleast 2 jets ==> events pass loose
    return veto = false;
  }
  return veto=true;
}


bool etau_analyzer::passDiElectronVeto(int eleIndex)
{
  std::vector<int> tmpCand; int tmpEleIndex1=-1; int tmpEleIndex2=-1;
  tmpCand.clear();
  bool veto = true;
  bool awayFromEverything=true;
  for(int iEle=0;iEle<nEle;iEle++)
    {
      bool kinematic = false;
      if( elePt->at(iEle) > 15
	  && fabs(eleEta->at(iEle))< 2.5
	  && fabs(eleD0->at(iEle)) < 0.045
	  && fabs(eleDz->at(iEle)) < 0.2
       	  ) kinematic = true;
      bool electronId =false;
      if( eleIDbit->at(iEle)>>3&1==1) electronId =true; // cut based electron id veto
      bool relative_iso = false;    
      float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
      if( relEleIso < 0.3 ) relative_iso = true;
      if( kinematic && electronId && relative_iso ){
	tmpCand.push_back(iEle);
      }	
    }
  std::vector<int> iElePlus; iElePlus.clear(); 
  std::vector<int> iEleMinus; iEleMinus.clear();
  for(int i=0; i<tmpCand.size(); i++){
    if(eleCharge->at(tmpCand[i]) < 0) iEleMinus.push_back(tmpCand[i]);
    if(eleCharge->at(tmpCand[i]) > 0) iElePlus.push_back(tmpCand[i]);
  }
  if(iElePlus.size()>0 && iEleMinus.size()>0){
    double deltaR= delta_R(elePhi->at(iEleMinus[0]), eleEta->at(iEleMinus[0]), elePhi->at(iElePlus[0]), eleEta->at(iElePlus[0]));
    if (deltaR > 0.15 && eleCharge->at(iElePlus[0])*eleCharge->at(iEleMinus[0])<0) {
      return false;
    }
  }
  return true;
  
}

bool etau_analyzer::eVetoZTTp001dxyz(int eleIndex, int tauIndex){
  std::vector<int> tmpCand;  tmpCand.clear();
  std::vector<int> output;  output.clear();
  bool awayFromEverything = true;   int tmpEleIndex=-1;
  //Loop over electrons      
  for(int iEle=0;iEle<nEle;iEle++)
    {
      if(iEle==eleIndex)continue;
      bool kinematic = false;
      if( elePt->at(iEle) > 10
	  && fabs(eleEta->at(iEle))< 2.5
	  && fabs(eleD0->at(iEle)) < 0.045
	  && fabs(eleDz->at(iEle)) < 0.2
	  && eleConvVeto->at(iEle)==1 && eleConvVeto->at(iEle)==1
       	  ) kinematic = true;
      bool electronId =false;
      if( eleIDbit->at(iEle)>>8&1==1) electronId =true; // cut based electron id veto
      bool relative_iso = false;    
      float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
      if( relEleIso < 0.3 ) relative_iso = true;
      if( kinematic && electronId && relative_iso ){
	tmpCand.push_back(iEle);
      }	
    }
  if(tmpCand.size()>0)
    {
      for(int i=0;i<tmpCand.size();i++)
	{
	  double deltaR_et = delta_R(tau_Phi->at(tauIndex), tau_Eta->at(tauIndex), elePhi->at(tmpCand[i]), eleEta->at(tmpCand[i]));
	  double deltaR_ee = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), elePhi->at(tmpCand[i]), eleEta->at(tmpCand[i]));
	  if(deltaR_et<0.5 && deltaR_ee<0.5)
	    output.push_back(i);
	}
    }
  if(output.size() >1 )
    return false;
  else
    return true;
    
   
}
bool etau_analyzer::mVetoZTTp001dxyz(int eleIndex, int tauIndex){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool awayFromEverything = true;   int tmpMuIndex=-1;
  //Loop over muons
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
  double deltaRm1=0.0; double deltaRm2=0.0;
  if(tmpCand.size() > 0 )
    { 
      deltaRm1 = delta_R(elePhi->at(eleIndex),eleEta->at(eleIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      deltaRm2 = delta_R(tau_Phi->at(tauIndex),tau_Eta->at(tauIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      if(deltaRm1>0.5 && deltaRm2>0.5 ){
	return false;
      }
    }
  else
    return true;

   
}
int etau_analyzer::myGenMaching(int tauIndex)
{
  double recotau_eta=tau_Eta->at(tauIndex);
  double recotau_phi=tau_Phi->at(tauIndex);
  double closestEle=999;  double closestMu=999;
  double closestETau=999;  double closestMTau=999;  double closestHTau=999;  double closestDR=999;
  double genLeptonEta=0;
  double genLeptonPhi=0;
  bool prompt_ele=false;  bool tau_ele=false; bool tau_mu=false; bool tau_tauh=false;
  bool prompt_mu=false;
  for(int imc=0; imc<nMC; imc++){
    genLeptonEta=mcEta->at(imc);
    genLeptonPhi=mcPhi->at(imc);
    double mc_tau_dr= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
    if(mc_tau_dr<closestDR)
      closestDR=mc_tau_dr;
  }


  for(int imc=0; imc<nMC; imc++){
    genLeptonEta=mcEta->at(imc);
    genLeptonPhi=mcPhi->at(imc);
    double dr_tau_lepton=dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
    //prompt_ele=false; prompt_mu=false; tau_ele=false; tau_mu=false; tau_tauh=false;
    
    ///// prompt electrons
    if(mcPt->at(imc)>8 && abs(mcPID->at(imc))==11 && mcStatusFlag->at(imc)>>1&1==1)
      {
	//dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
	if( dr_tau_lepton<0.2 && closestEle>dr_tau_lepton)
	  {closestEle=dr_tau_lepton; prompt_ele=true; }
	
      }
    ///// prompt muons
    if(mcPt->at(imc)>8 && abs(mcPID->at(imc))==13 && mcStatusFlag->at(imc)>>1&1==1)
      {
	//dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
        if( dr_tau_lepton<0.2 && closestMu>dr_tau_lepton)
          {closestMu=dr_tau_lepton; prompt_mu=true; }
      }
    ///// tau -> electrons
    if(mcPt->at(imc)>8 && abs(mcPID->at(imc))==11 && mcStatusFlag->at(imc)>>5&1==1)
      {
	//dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
        if( dr_tau_lepton<0.2 && closestETau>dr_tau_lepton)
          {closestETau=dr_tau_lepton;  tau_ele=true; }
      }
    ///// tau -> muons
    if(mcPt->at(imc)>8 && abs(mcPID->at(imc))==13 && mcStatusFlag->at(imc)>>5&1==1)
      {
	//dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
        if( dr_tau_lepton<0.2 && closestMTau>dr_tau_lepton)
          {closestMTau=dr_tau_lepton;  tau_mu=true; }
      }
    ///// tau -> tau hadronic
    if(mcPt->at(imc)>15 &&  abs(mcPID->at(imc))!=13 &&  abs(mcPID->at(imc))!=11 )
      {
	dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
        if( dr_tau_lepton<0.2 && closestHTau>dr_tau_lepton)
          {closestHTau=dr_tau_lepton;   tau_tauh=true; }
      } 
  }
  double closestLTau =  min(closestETau, closestMTau);
  if(closestHTau < closestLTau)
    closestLTau=closestHTau;
  //closestDR = min(closestLTau, min(closestEle, closestMu) );
  int genMatch=0;
  //cout<<"closestDR: "<<closestDR<<" closestEle:"<<closestEle<<" closestMu:"<<closestMu<<" closestETau:"<<closestETau<<" closestMTau:"<<closestMTau<<" closestHTau:"<<closestHTau<<endl;
  
  if( (prompt_ele || prompt_mu))
    {
      if(closestEle<0.2 && prompt_ele)
	//return 1;
	genMatch=1;
      else if(closestMu<0.2 && prompt_mu)				
	//return 2;
	genMatch=2;
    }
  else if(closestDR <= closestLTau)
    {
      if(closestETau<0.2 && closestETau< min(closestMTau, closestHTau) && tau_ele) //return 3;
	genMatch=3;
      else if(closestMTau<0.2 && closestMTau< min(closestETau, closestHTau) && tau_mu) //return 4;
	genMatch=4;
      else if(closestHTau<0.2 && closestHTau< min(closestETau, closestMTau) && tau_tauh) //return 5;
	genMatch=5;
    }
  else
    genMatch=6;


  return genMatch;

}
int etau_analyzer::myGenMaching1(int eleIndex)
{
  double recotau_eta=eleEta->at(eleIndex);
  double recotau_phi=elePhi->at(eleIndex);
  double closestEle=999;  double closestMu=999;
  double closestETau=999;  double closestMTau=999;  double closestHTau=999;  double closestDR=999;
  double genLeptonEta=0;
  double genLeptonPhi=0;
  bool prompt_ele=false;  bool tau_ele=false; bool tau_mu=false; bool tau_tauh=false;
  bool prompt_mu=false;
  for(int imc=0; imc<nMC; imc++){
    genLeptonEta=mcEta->at(imc);
    genLeptonPhi=mcPhi->at(imc);
    double mc_tau_dr= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
    if(mc_tau_dr<closestDR)
      closestDR=mc_tau_dr;
  }


  for(int imc=0; imc<nMC; imc++){
    genLeptonEta=mcEta->at(imc);
    genLeptonPhi=mcPhi->at(imc);
    double dr_tau_lepton=dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
    //prompt_ele=false; prompt_mu=false; tau_ele=false; tau_mu=false; tau_tauh=false;
    
    ///// prompt electrons
    if(mcPt->at(imc)>8 && abs(mcPID->at(imc))==11 && mcStatusFlag->at(imc)>>1&1==1)
      {
	//dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
	if( dr_tau_lepton<0.2 && closestEle>dr_tau_lepton)
	  {closestEle=dr_tau_lepton; prompt_ele=true; }
	
      }
    ///// prompt muons
    if(mcPt->at(imc)>8 && abs(mcPID->at(imc))==13 && mcStatusFlag->at(imc)>>1&1==1)
      {
	//dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
        if( dr_tau_lepton<0.2 && closestMu>dr_tau_lepton)
          {closestMu=dr_tau_lepton; prompt_mu=true; }
      }
    ///// tau -> electrons
    if(mcPt->at(imc)>8 && abs(mcPID->at(imc))==11 && mcStatusFlag->at(imc)>>5&1==1)
      {
	//dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
        if( dr_tau_lepton<0.2 && closestETau>dr_tau_lepton)
          {closestETau=dr_tau_lepton;  tau_ele=true; }
      }
    ///// tau -> muons
    if(mcPt->at(imc)>8 && abs(mcPID->at(imc))==13 && mcStatusFlag->at(imc)>>5&1==1)
      {
	//dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
        if( dr_tau_lepton<0.2 && closestMTau>dr_tau_lepton)
          {closestMTau=dr_tau_lepton;  tau_mu=true; }
      }
    ///// tau -> tau hadronic
    if(mcPt->at(imc)>15 &&  abs(mcPID->at(imc))!=13 &&  abs(mcPID->at(imc))!=11 )
      {
	dr_tau_lepton= dR(recotau_eta, recotau_phi, genLeptonEta, genLeptonPhi);
        if( dr_tau_lepton<0.2 && closestHTau>dr_tau_lepton)
          {closestHTau=dr_tau_lepton;   tau_tauh=true; }
      } 
  }
  double closestLTau =  min(closestETau, closestMTau);
  if(closestHTau < closestLTau)
    closestLTau=closestHTau;
  //closestDR = min(closestLTau, min(closestEle, closestMu) );
  int genMatch=0;
  //cout<<"closestDR: "<<closestDR<<" closestEle:"<<closestEle<<" closestMu:"<<closestMu<<" closestETau:"<<closestETau<<" closestMTau:"<<closestMTau<<" closestHTau:"<<closestHTau<<endl;
  
  if( (prompt_ele || prompt_mu))
    {
      if(closestEle<0.2 && prompt_ele)
	//return 1;
	genMatch=1;
      else if(closestMu<0.2 && prompt_mu)				
	//return 2;
	genMatch=2;
    }
  else if(closestDR <= closestLTau)
    {
      if(closestETau<0.2 && closestETau< min(closestMTau, closestHTau) && tau_ele) //return 3;
	genMatch=3;
      else if(closestMTau<0.2 && closestMTau< min(closestETau, closestHTau) && tau_mu) //return 4;
	genMatch=4;
      else if(closestHTau<0.2 && closestHTau< min(closestETau, closestMTau) && tau_tauh) //return 5;
	genMatch=5;
    }
  else
    genMatch=6;

  
  return genMatch;

}
std::vector<int> etau_analyzer::gen_matching(){
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
  
  // if(tmpGenMatch.size() >0 )
  //   tmpCand=tmpGenMatch[0];
  // return tmpCand; 
  return tmpGenMatch;
}
bool  etau_analyzer::found_GenMatch(int genTau)
{
  std::vector<int> v = gen_matching();
  if (std::find(v.begin(), v.end(), genTau) != v.end())
    return true;
  
  return false;
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
bool etau_analyzer::hasGenTau(){
  bool found_genTau=false;
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>5&1==1) {  found_genTau=true;}
  }
  return found_genTau;
}
float etau_analyzer::exponential(float x,float a,float b,float c) {
  return a * TMath::Exp(-b * x) + c;
}
double etau_analyzer::getScaleFactors(  double elept, double taupt, double eleeta, double taueta, int taudm, int tauGenMatch, bool isFakebkg)
{
  bool debug=false;
  double rv_sf=1.0;
  double eleRecoSF_corr=1.0;
  double eleEffSF_corr=1.0;
  double eletrgsf_tmp=1.0;
  double eletrgsf=1.0;
  double sf_tauTrg = 1.0; double sf_htt_workspace=1.0;
  double sf_Zvtx=1.0;
  double sf_tauidSF_m = 1.0;
  double sf_tauidSF_vvvl = 1.0;
  double sf_tauesSF = 1.0;
  double sf_fakeEle = 1.0; double sf_fakeMu = 1.0;
  double sf_fakeEleES = 1.0; double sf_fakeMuES = 1.0;
  double sf_taufesSF = 1.0;
  
  eleRecoSF_corr=h_eleRecoSF_highpt->GetBinContent(h_eleRecoSF_highpt->GetXaxis()->FindBin(eleeta),h_eleRecoSF_highpt->GetYaxis()->FindBin(elept));
  if (debug==true ) std::cout<<"eleRecoSF_corr =  "<< eleRecoSF_corr<<std::endl;
  eleEffSF_corr=h_eleIDSF->GetBinContent(h_eleIDSF->GetXaxis()->FindBin(eleeta),h_eleIDSF->GetYaxis()->FindBin(elept));
  if (debug==true ) std::cout<<"eleEffSF_corr =  "<< eleEffSF_corr<<std::endl;
  if (debug==true ) std::cout<<"This works line 269 "<<std::endl;
  

  if(  tauGenMatch>=5 )
    {
      sf_tauidSF_m = h_tauidSF_m->GetBinContent(h_tauidSF_m->GetXaxis()->FindBin(taudm));
      //sf_tauidSF_m = fn_tauIDSF_m->Eval(tau_Pt->at(reco_tau));
      sf_tauidSF_vvvl = h_tauidSF_vvvl->GetBinContent(h_tauidSF_vvvl->GetXaxis()->FindBin(taudm));
      //sf_tauidSF_vvvl = fn_tauIDSF_vvl->Eval(tau_Pt->at(reco_tau));
      sf_tauesSF = h_tauesSF->GetBinContent(h_tauesSF->GetXaxis()->FindBin(taudm));
    }
  
  if(tauGenMatch < 5 ) /// electrons to pass Deep
    {
      if(taudm==0)
  	{
  	  if(abs(taueta) < 1.479 ) sf_fakeEle=0.98;
  	  if(abs(taueta) > 1.479 ) sf_fakeEle=0.80;
  	}
      if(taudm==1)
  	{
  	  if(abs(taueta) < 1.479 ) sf_fakeEle=1.07;
  	  if(abs(taueta) > 1.479 ) sf_fakeEle=0.64;
  	}
    }
  
  double tauPtCheck=taupt;
  if(taupt > 450 ) tauPtCheck = 450;
  else if ( taupt < 20 )  tauPtCheck = 20;
  
  if(taudm==0)  sf_tauTrg= h_tauTrgSF_dm0->GetBinContent(h_tauTrgSF_dm0->GetXaxis()->FindBin(tauPtCheck));
  if(taudm==1)  sf_tauTrg= h_tauTrgSF_dm1->GetBinContent(h_tauTrgSF_dm1->GetXaxis()->FindBin(tauPtCheck));
  if(taudm==10) sf_tauTrg= h_tauTrgSF_dm10->GetBinContent(h_tauTrgSF_dm10->GetXaxis()->FindBin(tauPtCheck));
  if(taudm==11) sf_tauTrg= h_tauTrgSF_dm11->GetBinContent(h_tauTrgSF_dm11->GetXaxis()->FindBin(tauPtCheck));
  
  w->var("e_pt")->setVal(elept);
  w->var("e_eta")->setVal(eleeta);
  w->var("t_pt")->setVal(taupt);
  w->var("t_mvadm")->setVal(taudm);
  double e_trk_sf=w->function("e_trk_ratio")->getVal();
  double e_idiso_sf=w->function("e_idiso_ic_ratio")->getVal();
  double e_trg_sf=w->function("e_trg_ic_ratio")->getVal();
  double e_trg24_sf=w->function("e_trg_24_ic_ratio")->getVal();
  double t_trg_sf=w->function("t_trg_ic_deeptau_medium_mvadm_etau_ratio")->getVal();
  double t_deepid_tightvsele_sf=w->function("t_deeptauid_mvadm_medium_tightvsele")->getVal();
  double zptmass_weight = 1.0;
  int genZCand=getZCand();
  if(genZCand>-1){
    w->var("z_gen_pt")->setVal(mcPt->at(genZCand));
    w->var("z_gen_mass")->setVal(mcMass->at(genZCand));
    // cout<< "zptmass_weight_nom = "<<w->function("zptmass_weight_nom")->getVal();
    zptmass_weight=w->function("zptmass_weight_nom")->getVal();
    //cout<<"z_gen_pt="<<mcPt->at(genZCand)<<" z_gen_mass="<<mcMass->at(genZCand)<< "  zptmass_weight="<<zptmass_weight<<endl;
  }
  
  if( (HLTEleMuX>>3&1==1 || HLTEleMuX>>61&1==1 || HLTEleMuX>>5&1==1) && elept>28.0 )
    e_trg_sf=e_trg_sf;
  else
    e_trg_sf=1.0;
  if((HLTTau>>1&1 ==1 ) && elept<28.0 && taupt > 35.0)
    e_trg24_sf=e_trg24_sf;
  else
    e_trg24_sf=1.0; 

  sf_htt_workspace=  e_trk_sf * e_idiso_sf *  e_trg24_sf * e_trg_sf;
  //cout<<"sf_htt_workspace : "<<sf_htt_workspace<<endl;

  rv_sf = eleRecoSF_corr * eleEffSF_corr * sf_tauidSF_m *sf_tauTrg* sf_htt_workspace;
  if(isFakebkg)
    rv_sf=rv_sf*sf_tauidSF_vvvl;
  if(rv_sf>0)
    return rv_sf;
  else
    return 1.0;

}
bool etau_analyzer::MatchTriggerFilter(int eleIndex, int tauIndex)
{
  bool filterele24tau30=false;
  bool filterele27=false;
  bool filterele32=false; bool filterele35=false;
  //HLT_Ele27_WPTight_Gsf_v
  if(HLTEleMuX>>3&1 == 1 && ( eleFiredSingleTrgs->at(eleIndex)>>12&1==1 || eleFiredSingleTrgs->at(eleIndex)>>13&1==1 || eleFiredSingleTrgs->at(eleIndex)>>14&1==1) ) filterele27=true;
  //HLT_Ele32_WPTight_Gsf_v
  if(HLTEleMuX>>61&1 == 1 && (eleFiredSingleTrgs->at(eleIndex)>>1&1==1 || eleFiredSingleTrgs->at(eleIndex)>>12&1==1 || eleFiredSingleTrgs->at(eleIndex)>>13&1==1 || eleFiredSingleTrgs->at(eleIndex)>>14&1==1  ))filterele32=true;
  //HLT_Ele35_WPTight_Gsf_v
  if(HLTEleMuX>>5&1 == 1  && (eleFiredSingleTrgs->at(eleIndex)>>12&1==1 || eleFiredSingleTrgs->at(eleIndex)>>13&1==1 || eleFiredSingleTrgs->at(eleIndex)>>14&1==1 )) filterele35=true;
  //HLT_Ele24_eta2p1_WPTight_
  if( HLTTau>>1&1 == 1 && (eleFiredSingleTrgs->at(eleIndex)>>13&1==1 ||  eleFiredSingleTrgs->at(eleIndex)>>14&1==1 ))
    filterele24tau30=true;
  // if(!is_MC)
  //   return true;
  // else if(is_MC)
  //   {
  //     if( filterele24tau30 || filterele27 || filterele32 || filterele35)
  //    return true;
  //     else
  //    return false;
  //   }
  
  return true;
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

void etau_analyzer::fillHist( string histNumber , int eleIndex, int tauIndex, bool isFakeBkg, float event_weight){
  string hNumber = histNumber;
  TLorentzVector eleP4; TLorentzVector tauP4;
  TLorentzVector initialMet; TLorentzVector initialMetTau; TLorentzVector corrMet;
  if(is_MC)
    initialMet=MetRecoilCorrections(eleIndex, tauIndex);
  else
    initialMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  tauP4.SetPtEtaPhiE(tau_Pt->at(tauIndex),tau_Eta->at(tauIndex),tau_Phi->at(tauIndex), tau_Energy->at(tauIndex));
  eleP4.SetPtEtaPhiE(elePt->at(eleIndex), eleEta->at(eleIndex), elePhi->at(eleIndex), eleE->at(eleIndex));
  initialMetTau=initialMet+tauP4;

  if(is_MC)
    applyESCorrections(eleP4,tauP4, eleIndex, tauIndex, isFakeBkg,
		       eleP4,tauP4);
  corrMet = initialMetTau-tauP4;
  if (!is_MC)
    corrMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  plotFill("elePt_"+hNumber,  eleP4.Pt() , 38 , 24 , 100,  event_weight);
  plotFill("eleEta_"+hNumber, eleEta->at(eleIndex), 25, -2.5, 2.5,  event_weight);
  plotFill("elePhi_"+hNumber, elePhi->at(eleIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("eleDz_"+hNumber,  eleDz->at(eleIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("eleD0_"+hNumber,  eleD0->at(eleIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("electronID_"+hNumber, eleIDbit->at(eleIndex)>>8&1, 4, -2, 2,  event_weight); // electronID
  float relEleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (eleP4.Pt());
  plotFill("relEleIso_"+hNumber, relEleIso, 15, 0, 0.3,  event_weight);
  plotFill("eleCharge_"+hNumber, eleCharge->at(eleIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tauP4.Pt() , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 25, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), tau_Phi->at(tauIndex), tau_Eta->at(tauIndex));
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);
  double deltaPhi = DeltaPhi(elePhi->at(eleIndex), tau_Phi->at(tauIndex));
  double deltaEta = fabs(eleEta->at(eleIndex) - tau_Eta->at(tauIndex));
  plotFill("deltaPhi_"+hNumber, deltaPhi , 30, -3.14, 3.14,  event_weight);
  plotFill("deltaEta_"+hNumber, deltaEta ,  25, -2.5, 2.5,  event_weight);

  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(eleIndex, tauIndex);
  plotFill("nJet_"+hNumber,  jetCand.size() , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, corrMet.Pt() , 20, 0, 100,  event_weight);
  plotFill("metPhi_"+hNumber, corrMet.Phi() , 30, -3.14, 3.14,  event_weight);
  double mT_eleMet = TMass_F( eleP4.Pt(), eleP4.Phi(), corrMet.Pt(),corrMet.Phi()  );
  plotFill("mT_eleMet_"+hNumber, mT_eleMet , 30, 0, 150,  event_weight);

  double visMass_mutau = VisMass_F(tauP4, eleP4);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  double HiggsPt = pTvecsum_F(tauP4, eleP4, corrMet );
  plotFill("higgsPt_"+hNumber,HiggsPt , 25, 0, 250,  event_weight);

  double tot_tr_mass = TotTMass_F(tauP4, eleP4, corrMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  int triggerBin1, triggerBin2, triggerBin3, triggerBin4;
  triggerBin1=triggerBin2=triggerBin3=triggerBin4=0;
  if( HLTEleMuX>>5&1 == 1 )  triggerBin4=4;
  if( HLTEleMuX>>61&1 == 1 ) triggerBin3=3;
  if( HLTEleMuX>>3&1 == 1 )  triggerBin2=2;
  if( HLTTau>>1&1 == 1 )     triggerBin1=1;
  if(triggerBin1>0)
    plotFill("trigger_"+hNumber, triggerBin1 , 5, 0, 5,  event_weight);
  if(triggerBin2>0)
    plotFill("trigger_"+hNumber, triggerBin2 , 5, 0, 5,  event_weight);
  if(triggerBin3>0)
    plotFill("trigger_"+hNumber, triggerBin3 , 5, 0, 5,  event_weight);
  if(triggerBin4>0)
    plotFill("trigger_"+hNumber, triggerBin4 , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(is_MC){
    if(myGenMaching(tauIndex)==1) genMatchBin=1;
    else if(myGenMaching(tauIndex)==2) genMatchBin=2;
    else if(myGenMaching(tauIndex)==3) genMatchBin=3;
    else if(myGenMaching(tauIndex)==4) genMatchBin=4;
    else if(myGenMaching(tauIndex)==5) genMatchBin=5;
    else if(myGenMaching(tauIndex)==6) genMatchBin=6;
  }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);
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
void etau_analyzer::makeTestPlot( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight){
  string hNumber = histNumber;
  std::vector<int> tmpCand; tmpCand.clear();
  for(int iEle=0;iEle<nEle;iEle++)
    {
      tmpCand.push_back(iEle);
    }
  plotFill("elePt_"+hNumber,  elePt->at(tmpCand[0]) , 38 , 24 , 100,  event_weight);
  //cout<<"     elePt_"<<hNumber<<" = "<< elePt->at(tmpCand[0])<<endl;
}


TLorentzVector etau_analyzer::MetRecoilCorrections(int eleIndex, int tauIndex){
  //// met recoil correction
  TLorentzVector mymet;
  TLorentzVector BosonP4, nuP4, nuP4tmp;
  TLorentzVector nu1P4, gentau1P4;
  TLorentzVector nu2P4, gentau2P4;
  TLorentzVector visGenP4;
  for(int i=0; i<nMC; i++)
    {
      if(mcPID->at(i)==23)
	BosonP4.SetPtEtaPhiE(mcPt->at(i), mcEta->at(i) , mcPhi->at(i) , mcE->at(i) );
    }
  //visGenP4=BosonP4;
  if(BosonP4.Pt()==0)
    {
      for(int i=0; i<nMC; i++)
	{
	  if(mcPID->at(i)==15)
	    gentau1P4.SetPtEtaPhiE(mcPt->at(i), mcEta->at(i) , mcPhi->at(i) , mcE->at(i) );
	  if(mcPID->at(i)==-15)
	    gentau2P4.SetPtEtaPhiE(mcPt->at(i), mcEta->at(i) , mcPhi->at(i) , mcE->at(i) );
	}
      BosonP4=gentau1P4+gentau2P4;
    }
  visGenP4=BosonP4;
  for(int i=0; i<nMC; i++)
    {
      if(abs(mcPID->at(i))==16 || abs(mcPID->at(i))==14 || abs(mcPID->at(i))==12)
	{
	  nuP4tmp.SetPtEtaPhiE(mcPt->at(i), mcEta->at(i) , mcPhi->at(i) , mcE->at(i) );
	  visGenP4=visGenP4-nuP4tmp;
	}
    }
  // apply recoil corrections on event-by-event basis (Type I PF MET)
  float pfmet=pfMET; float pfmetPhi=pfMETPhi;
  float pfmetcorr_ex=pfmet*cos(pfmetPhi); float pfmetcorr_ey=pfmet*sin(pfmetPhi);
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(eleIndex, tauIndex);
  //RecoilCorrector recoilPFMetCorrector("sf_files/HTT-utilities/RecoilCorrections/data/Type1_PFMET_2017.root"); // Type I PF MET 2017 (ReReco)
  //RecoilCorrector recoilPFMetCorrector;
  recoilPFMetCorrector.CorrectByMeanResolution(pfmet*cos(pfmetPhi), // uncorrected type I pf met px (float)
					       pfmet*sin(pfmetPhi), // uncorrected type I pf met py (float)
					       BosonP4.Px(), // generator Z/W/Higgs px (float)
					       BosonP4.Py(), // generator Z/W/Higgs py (float)
					       visGenP4.Px(), // generator visible Z/W/Higgs px (float)
					       visGenP4.Py(), // generator visible Z/W/Higgs py (float)
					       jetCand.size(),  // number of jets (hadronic jet multiplicity) (int)
					       pfmetcorr_ex, // corrected type I pf met px (float)
					       pfmetcorr_ey  // corrected type I pf met py (float)
					       );
  

  mymet.SetPxPyPzE(pfmetcorr_ex,pfmetcorr_ey,0,sqrt(pfmetcorr_ex*pfmetcorr_ex + pfmetcorr_ey*pfmetcorr_ey));
  return mymet;
}

void etau_analyzer::applyESCorrections(TLorentzVector eleP4, TLorentzVector tauP4, int eleIndex, int tauIndex, bool isFakebkg, TLorentzVector& eleP4Corr, TLorentzVector& tauP4Corr)
{
  
  eleP4Corr = eleP4*(eleCalibE->at(eleIndex)/eleP4.E());
  
  //if(!isFakebkg)
    {
      if (myGenMaching(tauIndex)>=5 && tau_DecayMode->at(tauIndex)==0) tauP4Corr=tauP4*1.007;
      else if (myGenMaching(tauIndex)>=5 && tau_DecayMode->at(tauIndex)==1) tauP4Corr=tauP4*0.998;
      else if (myGenMaching(tauIndex)>=5 && tau_DecayMode->at(tauIndex)==10) tauP4Corr=tauP4*1.001;
      if (  (myGenMaching(tauIndex)==1 || myGenMaching(tauIndex)==3) && tau_DecayMode->at(tauIndex)==0 ) tauP4Corr=tauP4*1.003;
      else if ( (myGenMaching(tauIndex)==1 || myGenMaching(tauIndex)==3) && tau_DecayMode->at(tauIndex)==1) tauP4Corr=tauP4*1.036;
    }
  double sf_taufesSF=1.0;
  if(isFakebkg)
    {
      if(tau_DecayMode->at(tauIndex)==0 && abs(tau_Eta->at(tauIndex))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(1);
      if(tau_DecayMode->at(tauIndex)==0 && abs(tau_Eta->at(tauIndex))>1.4 )  sf_taufesSF = h_taufesSF->Eval(3);
      if(tau_DecayMode->at(tauIndex)==1 && abs(tau_Eta->at(tauIndex))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(5);
      if(tau_DecayMode->at(tauIndex)==1 && abs(tau_Eta->at(tauIndex))>1.4 )  sf_taufesSF = h_taufesSF->Eval(7);
    }
  tauP4Corr=tauP4* sf_taufesSF;
}
