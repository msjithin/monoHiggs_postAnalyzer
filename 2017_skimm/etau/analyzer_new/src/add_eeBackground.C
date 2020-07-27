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


void etau_analyzer::Loop_ee(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
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
   
   fileName->cd();
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
       if(isMC)
       	 weight=weight*prefiringweight;
       if( isGoodVtx==false ) continue;
       //if( noisyJet2017()==true ) continue;
       if( found_DYjet_sample && hasGenTau())
	 eleTau_selector=true;
       else if( found_DYjet_sample && !hasGenTau() )
	 eleTau_selector=false;
       else if ( !found_DYjet_sample )
	 eleTau_selector=true;

       if( sample.Contains("ee_") )
	 {
	   //if(!found_GenMatch(5) && !found_GenMatch(6) )
	   eleEle_selector=true;
	   //else
	   //eleEle_selector=false;
	 }
       else
	 eleEle_selector=false;

       eleTau_selector=true;
       eleEle_selector=true;
       ////Trigger bit selection
       
       if( (HLTEleMuX>>3&1 == 1 )      //HLT_Ele27_WPTight_Gsf_v
	   || (HLTEleMuX>>61&1 == 1)  //HLT_Ele32_WPTight_Gsf_v
	   || (HLTEleMuX>>5&1 == 1)   //HLT_Ele35_WPTight_Gsf_v
	   )
       	 passSingleTriggerPaths=true;  //
       
       if( ( HLTTau>>1&1 == 1 ) )      //HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1
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
		       if( tauCand.size()>0  ) 
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
				   //if ( MatchTriggerFilter(reco_ele[0], reco_tau[0]) )
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
					 if( thirdLeptonVeto(reco_ele[0] , reco_tau[0], reco_ele2[0])  )
					   {
					     nPassedThirdLepVeto_dyll+=event_weight;
					     makeTestPlot("g_dyll", 0,0,0,event_weight);
					     if( passBjetVeto(reco_ele[0] , reco_tau[0], reco_ele2[0]) )
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
	   ////// fake background region - antiisolated  DY Z->ll begin
	   if(debug)cout<<"moving to fake bkg DY->ll "<<endl;
	   event_weight=weight;
	   eleCand.clear(); ele2Cand.clear();  tauCand.clear();
	   if(metFilters==0 )
	     {
	       if(isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	       nMETFiltersPassed_dyll_fr+=event_weight;
	       if(  passSingleTriggerPaths || passCrossTrigger  )
		 {
		   nSingleTrgPassed_dyll_fr+=event_weight;
		   if(debug)cout<<"trigger selected"<<endl;
		   eleCand = getEleCand(24,2.1);  ///// ele selected
		   if( eleCand.size() >0 ) 
		     { 
		       nGoodMuonPassed_dyll_fr+=event_weight;
		       if(debug)cout<<"this worked Line 514"<<endl;
		       std::vector<int> iElePlus; iElePlus.clear(); 
		       std::vector<int> iEleMinus; iEleMinus.clear();
		       for(int i=0; i<eleCand.size(); i++){
		       	 if(eleCharge->at(eleCand[i]) < 0) iEleMinus.push_back(eleCand[i]);
		       	 if(eleCharge->at(eleCand[i]) > 0) iElePlus.push_back(eleCand[i]);
		       }
		       if(iElePlus.size()>0 && iEleMinus.size()>0)
		       	 {
		       	   //if( eleCharge->at(eleCand[0])*eleCharge->at(ele2Cand[0]) <0 )
		       	   {
		       	     tauCand = getAISRTauCand(30,2.3);
		       	     if( tauCand.size()>0 ) 
		       	       {
		       		 nGoodTauPassed_dyll_fr+=event_weight;
		       		 if(debug)cout<<"fr tau selection passed"<<endl;
		       		 reco_ele.clear();reco_tau.clear(); reco_ele2.clear();
		       		 if( elePt->at(iEleMinus[0]) > elePt->at(iElePlus[0])   ) { reco_ele=iEleMinus; reco_ele2=iElePlus; }
		       		 else { reco_ele=iElePlus; reco_ele2=iEleMinus;}
		       		 reco_tau=tauCand;
				 
				   //if ( MatchTriggerFilter(reco_ele[0], reco_tau[0]) )
				     {
				       //if ( eleCharge->at(reco_ele[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
				       {
					 nGoodMuTauPassed_dyll_fr+=event_weight;
					 electronP4.SetPtEtaPhiE(elePt->at(reco_ele[0]), eleEta->at(reco_ele[0]), elePhi->at(reco_ele[0]), eleE->at(reco_ele[0]));
					 tauP4.SetPtEtaPhiE(tau_Pt->at(reco_tau[0]), tau_Eta->at(reco_tau[0]), tau_Phi->at(reco_tau[0]), tau_Energy->at(reco_tau[0]));
					 
					 event_weight = event_weight* getFR(reco_tau[0]);
					 
					 /////
					 if(debug)cout<<" fake bkg sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , true , isMC , debug ) <<endl;
					 if(isMC) event_weight = event_weight * getScaleFactors( reco_ele[0] , reco_tau[0] , true , isMC , debug );
					 /////
					 if( thirdLeptonVeto(reco_ele[0] , reco_tau[0], reco_ele2[0])  )
					   {
					     nPassedThirdLepVeto_dyll_fr+=event_weight;
					     
					     if( passBjetVeto(reco_ele[0] , reco_tau[0], reco_ele2[0]) )
					       {
						 nPassedBjetVeto_dyll_fr+=event_weight;
					     
						 double deltaR = delta_R(elePhi->at(reco_ele[0]),eleEta->at(reco_ele[0]), elePhi->at(reco_ele2[0]),eleEta->at(reco_ele2[0]));
						 if(deltaR > 0.5 )
						   {
						     nDeltaRPassed_dyll_fr+=event_weight;
						     if(debug)cout<<"this worked Line 425"<<endl;
						     if(debug)cout<<"evnt_weight = "<<event_weight<<endl;
						     fillHist_dyll("5_dyll_fr", reco_ele[0], reco_ele2[0], reco_tau[0], event_weight, isMC);
						     double mT_muMet = TMass_F((elePt->at(reco_ele[0])),(elePhi->at(reco_ele[0])),pfMET,pfMETPhi  );
						     if( mT_muMet < 50 ) 
						       {
							 fillHist_dyll("6_dyll_fr", reco_ele[0], reco_ele2[0], reco_tau[0], event_weight, isMC);
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

   // fileName->cd();
   // map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
   // map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
   // for (; iMap1 != jMap1; ++iMap1)
   //   nplot1(iMap1->first)->Write();
}
