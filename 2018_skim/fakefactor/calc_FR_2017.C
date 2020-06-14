/////calc_FR_2017.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom calc_FR_2017 analyze
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
#define calc_FR_2017_cxx
#include "calc_FR_2017.h"
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
//#include "makeHisto.h" 

using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{
  TStopwatch sw;
  sw.Start();
  
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
  
  calc_FR_2017 t(argv[1],argv[2], isMC);
  t.Loop(maxEvents,reportEvery, SampleName , isMC);
  sw.Stop();
  sw.Print();
  return 0;
}

void calc_FR_2017::Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
{

  
  if (fChain == 0) return;
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
  int nBackupTriggerEvents, nBTMediumEvents, nBTMediumHLTsinglePhoEvents, nEffPhoptden, nEffPhoptnum, nEffMETden, nEffMETnum;
  nBackupTriggerEvents = nBTMediumEvents = nBTMediumHLTsinglePhoEvents = nEffPhoptden = nEffPhoptnum = nEffMETden = nEffMETnum = 0;
  double nHLTPassed,n_eventWeight, nSingleTrgPassed, nGoodMuonPassed,nGoodMuonPassed_2, nElectronPtPassed, nGoodTauPassed, nTauPtPassed,numberOfEvents,nMETPassed, nDPhiPassed, nqcdden,nqcdnum, nMETFiltersPassed,nLeptonVetoPassed,nPassedBjetVeto,nNoisyCrystals,nDPhiJetMETPassed, nGoodMuMuPassed,nDeltaRPassed,nPassedThirdLepVeto, nPassedHiggsPtcut, nPassedVisibleMasscut, nPassedMETcut, nFinal_afterSelections, nGoodMuonPassed_qcd, nGoodTauPassed_qcd, nGoodMuTauPassed_qcd, nDeltaRPassed_qcd, nPassedThirdLepVeto_qcd,nPassedBjetVeto_qcd, nPassedHiggsPtcut_qcd, nPassedVisibleMasscut_qcd, nPassedMETcut_qcd, nFinal_afterSelections_qcd, nVLooseTau, nTightTau ;
  nHLTPassed = n_eventWeight = nSingleTrgPassed = nGoodMuonPassed=nGoodMuonPassed_2 = nElectronPtPassed = nGoodTauPassed = nTauPtPassed= numberOfEvents = nMETPassed = nDPhiPassed = nqcdden= nqcdnum=nMETFiltersPassed= nLeptonVetoPassed=nPassedBjetVeto=nNoisyCrystals=nDPhiJetMETPassed= nGoodMuMuPassed = nDeltaRPassed= nPassedThirdLepVeto=nPassedHiggsPtcut=nPassedVisibleMasscut=nPassedMETcut=nFinal_afterSelections=nGoodMuonPassed_qcd=nGoodTauPassed_qcd=nGoodMuTauPassed_qcd=nDeltaRPassed_qcd=nPassedThirdLepVeto_qcd=nPassedBjetVeto_qcd=nPassedHiggsPtcut_qcd=nPassedVisibleMasscut_qcd=nPassedMETcut_qcd=nFinal_afterSelections_qcd=nVLooseTau=nTightTau=0;
  
  std::vector<int> muCand_1;
  muCand_1.clear();
  std::vector<int> muCand_2;
  muCand_2.clear();
  
  std::vector<int> tauCand;
  tauCand.clear();
  
  std::vector<int> tauCand_t;
  tauCand_t.clear();
  TString sample = TString(SampleName);
  
   Float_t met_bins[16]={0.0, 15, 30, 45, 60, 90, 105, 120, 135, 150,175., 190., 250., 400., 600.0,800.0};
   TH1F* h_Events_level= new TH1F("Events_level","Events at each level of selection",15,0,30);
   TH1F* h_Cutflow= new TH1F("Cutflow","Events at each level of selection",6,0,12);
   Float_t pt_bins[8]= {20, 25, 30, 35, 40, 50, 60, 120};
   
   TH1F* h_tauPt_test= new TH1F("tauPt_inc_test","tauPt_inc_test",7,pt_bins ); h_tauPt_test->Sumw2();

   TH1F* h_tauPt_inc_t= new TH1F("tauPt_inc_tight","tauPt_inc_tight",7,pt_bins ); h_tauPt_inc_t->Sumw2();
   TH1F* h_tauPt_inc_vl= new TH1F("tauPt_inc_vl","tauPt_inc_vl",7,pt_bins ); h_tauPt_inc_vl->Sumw2();
   
   TH1F* h_tauPt_dm0_t= new TH1F("tauPt_dm0_tight","tauPt_dm0_tight",7,pt_bins ); h_tauPt_dm0_t->Sumw2();
   TH1F* h_tauPt_dm1_t= new TH1F("tauPt_dm1_tight","tauPt_dm1_tight",7,pt_bins ); h_tauPt_dm1_t->Sumw2();
   TH1F* h_tauPt_dm10_t= new TH1F("tauPt_dm_10_tight","tauPt_dm10_tight",7,pt_bins ); h_tauPt_dm10_t->Sumw2();  
   TH1F* h_tauPt_dm11_t= new TH1F("tauPt_dm_11_tight","tauPt_dm11_tight",7,pt_bins ); h_tauPt_dm11_t->Sumw2();

   TH1F* h_tauPt_dm0_vl= new TH1F("tauPt_dm0_veryloose","tauPt dm0  very loose",7,pt_bins ); h_tauPt_dm0_vl->Sumw2();
   TH1F* h_tauPt_dm1_vl= new TH1F("tauPt_dm1_veryloose","tauPt dm1 very loose",7,pt_bins ); h_tauPt_dm1_vl->Sumw2();
   TH1F* h_tauPt_dm10_vl= new TH1F("tauPt_dm10_veryloose","tauPt dm10 very loose",7,pt_bins ); h_tauPt_dm10_vl->Sumw2();
   TH1F* h_tauPt_dm11_vl= new TH1F("tauPt_dm11_veryloose","tauPt dm11 very loose",7,pt_bins ); h_tauPt_dm11_vl->Sumw2();

   TH1F* h_tauFF_dm0= new TH1F("tauFF_dm0","tau fake factor dm0",7,pt_bins ); h_tauFF_dm0->Sumw2();
   TH1F* h_tauFF_dm1= new TH1F("tauFF_dm1","tau fake factor dm1",7,pt_bins ); h_tauFF_dm1->Sumw2();
   TH1F* h_tauFF_dm10= new TH1F("tauFF_dm10","tau fake factor dm10",7,pt_bins ); h_tauFF_dm10->Sumw2();
   TH1F* h_tauFF_dm11= new TH1F("tauFF_dm11","tau fake factor dm11",7,pt_bins ); h_tauFF_dm11->Sumw2();



   TH1F* h_tauM_t = new TH1F("tauM_t","tau mass tight",20, 0.0, 2.0 ); h_tauM_t->Sumw2();
   TH1F* h_tauM_vl = new TH1F("tauM_vl","tau mass vl ",20, 0.0, 2.0 ); h_tauM_vl->Sumw2();

   TH1F* h_mumu_invM= new TH1F("mumu_invM","mu-mu invM",30, 0, 150 ); h_mumu_invM->Sumw2();
   TH1F* h_mutau_invM= new TH1F("mutau_invM","mu-tau invM",30, 0, 150 ); h_mutau_invM->Sumw2();
  
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
  
   //newtree=fChain->CloneTree(0);
   
   for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
     {
       nInspected+=1;
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       int report_i=0;
       double event_weight=1.0;
	  if(debug==true)std::cout<<"This works! line 154" <<endl;
	  int skip_i = -9;
	  // muCand_1 = getMuCand_1(30,2.1); // Muon pt, eta cuts
	  //if (muCand_1.size()>0) skip_i=muCand_1[0];
	  //muCand_2 = getMuCand_2(30,2.1, skip_i); // Muon pt, eta cuts
	  
	  //tauCand = getTauCand(20,2.3); // Tau pt, eta cuts
	  //particleCand = getParticleCand();
	  //	  if (run > 274421 && (event % 4) != 0) continue;
	  //fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;   
	  numberOfEvents+=event_weight;
	  if(debug==true)std::cout<<"This works! line 165" <<endl;

	  if(metFilters==0) 
	    {	    
	      if(debug==true)std::cout<<"This works! line 188" <<endl;
	      nMETFiltersPassed+=event_weight;
	      if(HLTEleMuX>>21&1 == 1 )//|| HLTEleMuX>>20&1 == 1 || HLTEleMuX>>32&1) // Single muon trigger: HLT_IsoMu27_v, HLT_IsoTkMu24_v , HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v
		{
		  nSingleTrgPassed+=event_weight;
		  tauCand_t = getTauCand(20,2.3);
		  
		  muCand_1 = getMuCand_1(30,2.4);
		  if(muCand_1.size() >0) 
		    {
		      if(debug==true)std::cout<<"This works! line 175" <<endl;
		      nGoodMuonPassed+=event_weight;
		      if (muCand_1.size()>0) skip_i=muCand_1[0];
		      muCand_2 = getMuCand_2(30,2.4, skip_i); // Muon pt, eta cuts
		      
		      if( muCand_2.size() >0 )
			{
			  nGoodMuonPassed_2 +=event_weight;
			  //std::cout<<"muon 1 pt = "<< muPt->at(muCand_1[0]) << " muon 2 pt = " << muPt->at(muCand_2[0]) <<endl;
			  if(debug==true)std::cout<<"This works! line 181" <<endl;
			  //std::cout<< "number of muons after selecting muons= "<< nMu <<endl;
			  if( (muCharge->at(muCand_1[0]))*(muCharge->at(muCand_2[0])) < 0 ) // opposite charge condition
			    {
			      if(debug==true)std::cout<<"This works! line 184" <<endl; nGoodMuMuPassed+=event_weight;
			      double deltaR_muMu = dR(muCand_1[0],muCand_2[0]);
			      if(deltaR_muMu > 0.5 )
				{ 
				  if(debug==true)std::cout<<"This works! line 188" <<endl;
				  nDeltaRPassed+=event_weight;
				  if( passBjetVeto() == true)
				    { nPassedBjetVeto+=event_weight;
				      if(debug==true)std::cout<<"This works! line 192" <<endl;
				      TLorentzVector myMu_1;
				      myMu_1.SetPtEtaPhiE(muPt->at(muCand_1[0]),muEta->at(muCand_1[0]),muPhi->at(muCand_1[0]), muE->at(muCand_1[0]));
				      TLorentzVector myMu_2;
				      myMu_2.SetPtEtaPhiE(muPt->at(muCand_2[0]),muEta->at(muCand_2[0]),muPhi->at(muCand_2[0]), muE->at(muCand_2[0]));
				      
				      double visMass_muMu = VisMass_F(myMu_1, myMu_2);
				      h_mumu_invM->Fill(visMass_muMu);
				      if( visMass_muMu > 80.0 && visMass_muMu< 100 ) 
					{
					  if(debug==true)std::cout<<"This works! line 201" <<endl;
					  nPassedVisibleMasscut+=event_weight;
					  //if (thirdLeptonVeto(muCand_1[0],muCand_2[0]) == true)
					  {
					    
					    tauCand = getTauCand(20,2.3); // Tau pt, eta cuts
					    
					    if( tauCand.size() >0 )
					      {
						nGoodTauPassed+=event_weight;
						TLorentzVector myTau;
						myTau.SetPtEtaPhiE(tau_Pt->at(tauCand[0]),tau_Eta->at(tauCand[0]),tau_Phi->at(tauCand[0]), tau_Energy->at(tauCand[0]));
						double visMass_muTau_1 = VisMass_F(myMu_1, myTau);
						double visMass_muTau_2 = VisMass_F(myMu_2, myTau);
						
						double deltaR_muTau_1 = dR_muTau(muCand_1[0],tauCand[0]);
						double deltaR_muTau_2 = dR_muTau(muCand_2[0],tauCand[0]);
						
						bool mu1_tau=false; bool mu2_tau=false;
						//if (!(deltaR_muTau_1 > 0.5 || deltaR_muTau_2 > 0.5 ))continue;
						if (muCharge->at(muCand_2[0])* tau_Charge->at(tauCand[0]) < 0
						    && !(deltaR_muTau_2 > 0.5))continue;
						if (muCharge->at(muCand_1[0])* tau_Charge->at(tauCand[0]) < 0
                                                    && !(deltaR_muTau_1 > 0.5))continue;
						
						//if (!(deltaR_muTau_1 > 0.5 || deltaR_muTau_2 > 0.5 ))continue;
						if (muCharge->at(muCand_2[0])* tau_Charge->at(tauCand[0]) < 0  
						     && visMass_muTau_2 > 80 && visMass_muTau_2<100 ) continue;
						if (muCharge->at(muCand_1[0])* tau_Charge->at(tauCand[0]) < 0
						    && visMass_muTau_1 > 80 && visMass_muTau_1<100 ) continue;
						
						
						if (tau_byVVVLooseDeepTau2017v2p1VSjet->at(tauCand[0])==1  )
						  {
						    nVLooseTau+=event_weight;
						    h_tauPt_inc_vl->Fill(tau_Pt->at(tauCand[0]), event_weight);
						    h_tauM_vl->Fill(myTau.M());
						    if (tau_DecayMode->at(tauCand[0])==0)
						      h_tauPt_dm0_vl->Fill(tau_Pt->at(tauCand[0]), event_weight);
						    if (tau_DecayMode->at(tauCand[0])==1)
						      h_tauPt_dm1_vl->Fill(tau_Pt->at(tauCand[0]), event_weight);
						    if (tau_DecayMode->at(tauCand[0])==10)
						      h_tauPt_dm10_vl->Fill(tau_Pt->at(tauCand[0]), event_weight);
						    if (tau_DecayMode->at(tauCand[0])==11)
                                                      h_tauPt_dm11_vl->Fill(tau_Pt->at(tauCand[0]), event_weight);

						  }
						
						if (tau_byMediumDeepTau2017v2p1VSjet->at(tauCand[0])==1  )
						  {
						    nTightTau+=event_weight;
						    h_tauPt_inc_t->Fill(tau_Pt->at(tauCand[0]), event_weight);
						    h_tauM_t->Fill(myTau.M());
						    if (tau_DecayMode->at(tauCand[0])==0)
						      h_tauPt_dm0_t->Fill(tau_Pt->at(tauCand[0]), event_weight);
						    if (tau_DecayMode->at(tauCand[0])==1)
						      h_tauPt_dm1_t->Fill(tau_Pt->at(tauCand[0]), event_weight);
						    if (tau_DecayMode->at(tauCand[0])==10)
						      h_tauPt_dm10_t->Fill(tau_Pt->at(tauCand[0]), event_weight);
						    if (tau_DecayMode->at(tauCand[0])==11)
                                                      h_tauPt_dm11_t->Fill(tau_Pt->at(tauCand[0]), event_weight);

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
   h_tauFF_dm0->Add(h_tauPt_dm0_t);
   h_tauFF_dm0->Divide(h_tauPt_dm0_vl);
   
   h_tauFF_dm1->Add(h_tauPt_dm1_t);
   h_tauFF_dm1->Divide(h_tauPt_dm1_vl);
   
   h_tauFF_dm10->Add(h_tauPt_dm10_t);
   h_tauFF_dm10->Divide(h_tauPt_dm10_vl);
   
   h_nEvents->SetBinContent(1, nInspected);
   cout<<"nVLooseTau = "<<nVLooseTau<<endl;
   cout<<"nTightTau  = "<<nTightTau<<endl;
}

void calc_FR_2017::BookHistos(const char* file2)
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


void calc_FR_2017::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{
  
  //h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);

}




std::vector<int> calc_FR_2017::getMuCand_1(double muPtCut, double muEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over muons                                                                     
  for(int iMu=0;iMu<nMu;iMu++)
    {
      bool kinematic = false;
      if( (*muPt)[iMu] > muPtCut  && fabs((*muEta)[iMu])< muEtaCut 
	  && fabs(muDz->at(iMu)) < 0.2 && fabs(muD0->at(iMu))<0.045 ) kinematic = true;
      bool muonId =false;
      if( muIDbit->at(iMu)>>1&1==1) muonId =true;
      bool relative_iso = false;
      //float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
      //if( relMuIso < 0.15 ) relative_iso = true;
      if( muIDbit->at(iMu)>>8&1==1 ) relative_iso = true; // PFMedium isolation
      if(muonId && kinematic && relative_iso){
	tmpCand.push_back(iMu);
      }                                                                                      
    }                                                                                       
  return tmpCand;
}
std::vector<int> calc_FR_2017::getMuCand_2(double muPtCut, double muEtaCut, int skip_i){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over muons
  for(int iMu=0;iMu<nMu;iMu++)
    {
      if (iMu==skip_i) continue;

      bool kinematic = false;
      if( (*muPt)[iMu] > muPtCut  && fabs((*muEta)[iMu])< muEtaCut 
	  && fabs(muDz->at(iMu)) < 0.2 && fabs(muD0->at(iMu))<0.045 ) kinematic = true;
      
      bool muonId =false;
      if( muIDbit->at(iMu)>>1&1==1) muonId =true;
      bool relative_iso = false;

      //float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
      //if( relMuIso < 0.15 ) relative_iso = true;
      if( muIDbit->at(iMu)>>8&1==1 ) relative_iso = true; // PFMedium isolation
      if(muonId && kinematic && relative_iso){
        tmpCand.push_back(iMu);
      }
    }
  return tmpCand;
}


std::vector<int> calc_FR_2017::getTauCand(double tauPtCut, double tauEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over taus      
  for(int iTau=0;iTau<nTau;iTau++)
    {
      bool kinematic = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool isdecayModeFindingNewDMs = false;
      if( tau_Pt->at(iTau) > tauPtCut  && fabs( tau_Eta->at(iTau))< tauEtaCut && fabs(tau_LeadChargedHadron_dz->at(iTau))<0.2 )kinematic = true;
      if( tau_byTightDeepTau2017v2p1VSmu->at(iTau) == 1 && tau_byVVLooseDeepTau2017v2p1VSe->at(iTau)==1) tauId = true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) isdecayModeFindingNewDMs = true;
      //if( tauDecayMode->at(iTau)==0 || tauDecayMode->at(iTau)==1 || tauDecayMode->at(iTau)==10 ) decayModeCut = true;
      
      if(kinematic == true && tauId==true && isdecayModeFindingNewDMs==true )
	//if(tauId==true && kinematic==true && decayModeCut==true)
	{
	  tmpCand.push_back(iTau);
	}                                                           
    }                                                                                       
  return tmpCand;
}



bool calc_FR_2017::thirdLeptonVeto(int muon_1, int muon_2)
{
  //  std::vector<int> tmpCand;
  //tmpCand.clear();
  int tmpCand = 0;
  bool thirdLep = false;
  bool thirdLepVeto=true;
   //Loop over muons                                                                                                                       
   for(int iMu=0; iMu < nMu;iMu++)
     {
       if (iMu==muon_1 || iMu==muon_2) continue;
       bool kinematic = false;
       if( (*muPt)[iMu] > 10.0  && fabs((*muEta)[iMu])< 2.4 ) kinematic = true;
       bool muonId =false;
       if( muIDbit->at(iMu)>>6&1==1) muonId =true;
       bool relative_iso = false;

       //float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
       //if( relMuIso < 0.3 ) relative_iso = true;
       if( muIDbit->at(iMu)>>16&1==1) relative_iso =true;

       if( kinematic==true && muonId==true && relative_iso == true){
	 //tmpCand.push_back(iMu);
	 tmpCand++;
       }                                                                                                                                                    
     }                                                                      
   if(tmpCand > 0){ thirdLep = true; thirdLepVeto=false;}
   //   if( thirdLep = true ) thirdLepVeto=false;
   //if( tmpCand = 0 ) thirdLepVeto=true;

   return thirdLepVeto;

}




//Veto failed if a electron is found that passes Loose Muon ID, Loose Muon Isolation, and elePtcut, and does not overlap the candidate photon within dR of 0.5                                                                                                                                                                    

double calc_FR_2017::dR(int mu_1_index, int mu_2_index )
{
  double deltaeta = abs(muEta->at(mu_1_index) - muEta->at(mu_2_index));
  double muon_1_Phi = muPhi->at(mu_1_index);
  double muon_2_Phi = muPhi->at(mu_2_index);

  double deltaphi = DeltaPhi(muon_1_Phi, muon_2_Phi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}



double calc_FR_2017::dR_muTau(int mu_2_index, int tau_index )
{
  double deltaeta = abs(muEta->at(mu_2_index) - tau_Eta->at(tau_index));
  double muon_2_Phi = muPhi->at(mu_2_index);
  double tauPhi = tau_Phi->at(tau_index);

  double deltaphi = DeltaPhi(muon_2_Phi, tauPhi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

double calc_FR_2017::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
}

float calc_FR_2017::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
    return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

float calc_FR_2017::TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float totalTMass = (a + b+ met).M();
  return totalTMass;
}


float calc_FR_2017::VisMass_F(TLorentzVector a, TLorentzVector b){
  float visibleMass = (a + b).M();
  return visibleMass;
}

float calc_FR_2017::pTvecsum_F(float pt1, float pt2, float phi1, float phi2) {
  float pt_vecSum = sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2), 2) + pow(pt1*sin(phi1) + pt2*sin(phi2), 2));
  return pt_vecSum;
}

bool calc_FR_2017::passBjetVeto()
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
