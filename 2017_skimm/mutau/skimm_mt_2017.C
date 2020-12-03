/////skimm_mt_2017.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom skimm_mt_2017 analyze
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
#define skimm_mt_2017_cxx
#include "skimm_mt_2017.h"
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

using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{
  TStopwatch sw;
  sw.Start();
  
  myMap1 = new map<string, TH1F*>();
  //myMap2 = new map<string, TH2F*>();
  std::string SampleName = argv[5];
  std::string isMC  = argv[6];
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
  
  skimm_mt_2017 t(argv[1],argv[2], isMC);
  t.Loop(maxEvents,reportEvery, SampleName , isMC);
  delete myMap1;
  sw.Stop();
  sw.Print();
  return 0;
}

void skimm_mt_2017::Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
{

  
  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;
  int report_=0;
  int report_test=0;
  int nInspected;
  nInspected = 0;
  float nPassedSkimmed=0;
  float triggerPassed=0;
  float trgFilterPass=0;
  double nInspected_genWeighted;  
  nInspected_genWeighted = 0.0; 
  bool debug=false;  double netWeight = 1.0;
   TString sample = TString(SampleName);
  
   bool fill_hist = false;
   bool isMC = false;
   if( _isMC_=="MC" ) { isMC=true; fill_hist=true; }
   else if ( _isMC_=="data" ) { isMC=false; fill_hist=false; }
   //bool debug=true;

   Long64_t nentries = fChain->GetEntries();
   if ( isMC==true ) std::cout<<".... MC file ..... "<<std::endl;
   std::cout<<"Coming in: "<<std::endl;
   std::cout<<"nentries:"<<nentries<<std::endl;
   //Look at up to maxEvents events, or all if maxEvents == -1.
   Long64_t nentriesToCheck = nentries;
   if (maxEvents != -1LL && nentries > maxEvents)
     nentriesToCheck = maxEvents;
   nTotal = nentriesToCheck;
   Long64_t nbytes = 0, nb = 0;
   
   std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  
   newtree=fChain->CloneTree(0);
   
   for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
     {
       nInspected+=1;
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
        int report_i=0;
       
       if(skimming_Htt()==true)
	 {
	   nPassedSkimmed+=1;
	   //fillOutTree();
	   newtree->Fill();
	 }
       // for(int itrg=0; itrg<13; itrg++)
       // 	 {
       // 	   if(HLTTau>>itrg&1==1 )
       // 	     {
       // 	       cout<<" "<<endl;
       // 	       cout<<"HLTTau "<<itrg<<"  event:"<<jentry<<endl;
       // 	       for(int imu=0; imu<nMu; imu++)
       // 		 {
       // 		   cout<<" "<<endl;
       // 		   cout<<"  muon:"<<imu<<endl;
       // 		   for(int i=0; i<56; i++)
       // 		     {
       // 		       if(muFiredTrgs->at(imu)>>i&1==1)
       // 			 cout<<"  filter: "<<i<<" ";
       // 		     }
       // 		 }
       // 	       for(int itau=0; itau<nTau; itau++)
       //           {
       //             cout<<" "<<endl;
       //             cout<<"  tau:"<<itau<<endl;
       //             for(int i=0; i<18; i++)
       //               {
       //                 if(tauFiredTrgs->at(itau)>>i&1==1)
       //                   cout<<"  filter: "<<i<<" ";
       //               }
       //           }
       // 	     }
       // 	 }
       //if(HLTTau>>0&1==1 )
       // if(HLTEleMuX>>60&1==1)
       // 	 {
       // 	   //cout<<" "<<endl;
       // 	   triggerPassed++;
       // 	   for(int imu=0; imu<nMu; imu++)
       // 	     {
       // 	       //cout<<" "<<endl;
       // 	       //cout<<"  muon:"<<imu<<endl;
       // 	       for(int i=33; i<56; i++)
       // 		 {
       // 		   if(muFiredTrgs->at(imu)>>i&1==1)
       // 		     {
       // 		       trgFilterPass++;
       // 		       //cout<<"  filter: "<<i<<" ";
       // 		     }
       // 		 }
       // 	     }
       // 	 }
       
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
   h_nEvents->SetBinContent(1, nInspected);
   cout<<"nPassedSkimmed = "<<nPassedSkimmed<<endl;
   cout<<"triggerPassed  = "<<triggerPassed<<endl;
   cout<<"trgFilterPass  = "<<trgFilterPass<<endl;
   cout<<"triggerPassed*23  = "<<triggerPassed*23<<endl;
}

void skimm_mt_2017::BookHistos(const char* file2)
{
  fileName = new TFile(file2, "RECREATE");
  h_nEvents=new TH1F("nEvents", "nEvents", 5, 0, 5);
  //tree = new TTree("eventTree","eventTree");
  //tree->Branch("run_",&run_);
  //tree->Branch("event_",&event_); 
  //tree->Branch("lumis_",&lumis_); 
  //makeOutputTree(tree);
  fileName->cd();
  
}

//Fill the sequential histos at a particular spot in the sequence


void skimm_mt_2017::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{
  
  //h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);

}




//---------------------------------------------------                                                                                                                                
// get a electron candiate based on pt eta and isolation                                                                                                                               
//----------------------------------------------------                                                                                                                               

bool skimm_mt_2017::skimming_Htt(){
  bool tmpCand=false;
  bool muFound=false; bool tauFound=false; bool drCutPassed=false;
  std::vector<int> tmpMuCand;    tmpMuCand.clear();
  std::vector<int> tmpTauCand;    tmpTauCand.clear();

  for(int iMu=0; iMu<nMu;iMu++){
    float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
    
    if(fabs(muDz->at(iMu)) < 0.2 && 
       fabs(muD0->at(iMu))<0.045 && 
       muPt->at(iMu) > 19.5      &&
       fabs(muEta->at(iMu))< 2.4 &&
       ( relMuIso < 0.3 ||  muIDbit->at(iMu)>>0&1==1 || muIDbit->at(iMu)>>6&1==1 || muIDbit->at(iMu)>>16&1==1  )
       ) {
      muFound=true;
      tmpMuCand.push_back(iMu);
    }
  }
  
  for(int iTau=0; iTau<nTau;iTau++){
    if( tau_Pt->at(iTau) > 19.5           
	&& fabs( tau_Eta->at(iTau))< 2.3     
	&& (tau_DecayMode->at(iTau) !=5 && tau_DecayMode->at(iTau)!=6) 
	&& (tau_IDbits->at(iTau)>>2&1==1 || tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1 ) 
	&& (tau_IDbits->at(iTau)>>4&1==1 || tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 ) 
	&& (tau_IDbits->at(iTau)>>13&1==1 || tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 )
	) {
      tauFound = true;
      tmpTauCand.push_back(iTau);
    }
  }
  
  for(int iMu=0; iMu<tmpMuCand.size();iMu++){
    for(int iTau=0; iTau<tmpTauCand.size();iTau++){
      double deltaR = dR(tmpMuCand[iMu], tmpTauCand[iTau]); 
      if(deltaR > 0.3 )
	drCutPassed=true;
    }
  }
  if(muFound==true && tauFound == true && drCutPassed==true)
    tmpCand=true;
  return tmpCand;
}

std::vector<int> skimm_mt_2017::skimmed_Mu(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int iMu=0; iMu<nMu;iMu++){
    float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
    if(fabs(muDz->at(iMu)) < 0.2 &&
       fabs(muD0->at(iMu))<0.045 &&
       muPt->at(iMu) > 19.5      &&
       fabs(muEta->at(iMu))< 2.4 &&
       ( muIDbit->at(iMu)>>0&1==1  || relMuIso < 0.30)
       ) {
      tmpCand.push_back(iMu);
    }
  }
  return tmpCand; 
}

std::vector<int> skimm_mt_2017::skimmed_Tau(){
  std::vector<int> tmpCand;    tmpCand.clear();
  for(int iTau=0; iTau<nTau;iTau++){
    if( tau_Pt->at(iTau) > 19.5           &&
	fabs( tau_Eta->at(iTau))< 2.3     &&
	(tau_DecayMode->at(iTau) !=5 && tau_DecayMode->at(iTau)!=6) &&
	(tau_IDbits->at(iTau)>>2&1==1 || tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1 ) &&
	(tau_IDbits->at(iTau)>>4&1==1 || tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 )&&
	(tau_IDbits->at(iTau)>>13&1==1 || tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1  )
	) {
      tmpCand.push_back(iTau);
    }
  }
  return tmpCand;
}

void skimm_mt_2017::fillOutTree(){
  run_ = run;
  event_= event;
  lumis_=lumis;
  tree->Fill();
}



double skimm_mt_2017::dR(int mu_index, int tau_index)
{
  double deltaeta = abs(muEta->at(mu_index) - tau_Eta->at(tau_index));
  double muonPhi = muPhi->at(mu_index);
  double tauPhi = tau_Phi->at(tau_index);

  double deltaphi = DeltaPhi(muonPhi, tauPhi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}

double skimm_mt_2017::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
}
