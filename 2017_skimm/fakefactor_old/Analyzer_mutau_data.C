////Analyzer_mutau_data.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom Analyzer_mutau_data analyze
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
#define Analyzer_mutau_data_cxx
#include "Analyzer_mutau_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
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
#include <set>
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
#include <TLorentzVector.h>

using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{

  std::string SampleName = argv[7];

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
  Analyzer_mutau_data t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery, SampleName);
  std::cout<<"***************************"<<std::endl;
  std::cout<<"Output written to  "<<argv[2]<<std::endl;
  return 0;
}

void Analyzer_mutau_data::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
{

   if (fChain == 0) return;
   int nTotal;
   nTotal = 0;

   int nInspected;
   nInspected = 0;
   double nInspected_genWeighted;  
   nInspected_genWeighted = 0.0; 
   bool debug=false;  double netWeight = 1.0;
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

   int iphi = 41;
   int ieta = 5;


   Long64_t nentries = fChain->GetEntries();
   std::cout<<"Coming in: "<<std::endl;
   std::cout<<"nentries:"<<nentries<<std::endl;
   //Look at up to maxEvents events, or all if maxEvents == -1.
   Long64_t nentriesToCheck = nentries;
	if (maxEvents != -1LL && nentries > maxEvents)
		nentriesToCheck = maxEvents;
	nTotal = nentriesToCheck;
	Long64_t nbytes = 0, nb = 0;

	std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
	TStopwatch sw;
	sw.Start();
	for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
	{
	  
	  event_.clear();
	  event_info.clear();
	  
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;
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
	      nMETFiltersPassed+=event_weight;
	      if(HLTEleMuX>>19&1 == 1 )//|| HLTEleMuX>>20&1 == 1 || HLTEleMuX>>32&1) // Single muon trigger: HLT_IsoMu27_v, HLT_IsoTkMu24_v , HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v
		{
		  nSingleTrgPassed+=event_weight;
		  tauCand_t = getTauCand(20,2.3);
		  if (tauCand_t.size() >0)
		    if ( taubyVLooseIsolationMVArun2017v2DBoldDMwLT2017->at(tauCand_t[0])==1 ) 
		      h_tauPt_test->Fill(tauPt->at(tauCand_t[0]), event_weight);
		  
		  
		  //
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
				      myMu_1.SetPtEtaPhiE(muPt->at(muCand_1[0]),muEta->at(muCand_1[0]),muPhi->at(muCand_1[0]), muEn->at(muCand_1[0]));
				      TLorentzVector myMu_2;
				      myMu_2.SetPtEtaPhiE(muPt->at(muCand_2[0]),muEta->at(muCand_2[0]),muPhi->at(muCand_2[0]), muEn->at(muCand_2[0]));
				      
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
						myTau.SetPtEtaPhiE(tauPt->at(tauCand[0]),tauEta->at(tauCand[0]),tauPhi->at(tauCand[0]), tauEnergy->at(tauCand[0]));
						double visMass_muTau_1 = VisMass_F(myMu_1, myTau);
						double visMass_muTau_2 = VisMass_F(myMu_2, myTau);
						
						double deltaR_muTau_1 = dR_muTau(muCand_1[0],tauCand[0]);
						double deltaR_muTau_2 = dR_muTau(muCand_2[0],tauCand[0]);
						
						bool mu1_tau=false; bool mu2_tau=false;
						//if (!(deltaR_muTau_1 > 0.5 || deltaR_muTau_2 > 0.5 ))continue;
						if (muCharge->at(muCand_2[0])* tauCharge->at(tauCand[0]) < 0
						    && !(deltaR_muTau_2 > 0.5))continue;
						if (muCharge->at(muCand_1[0])* tauCharge->at(tauCand[0]) < 0
                                                    && !(deltaR_muTau_1 > 0.5))continue;
						
						//if (!(deltaR_muTau_1 > 0.5 || deltaR_muTau_2 > 0.5 ))continue;
						if (muCharge->at(muCand_2[0])* tauCharge->at(tauCand[0]) < 0  
						     && visMass_muTau_2 > 80 && visMass_muTau_2<100 ) continue;
						if (muCharge->at(muCand_1[0])* tauCharge->at(tauCand[0]) < 0
						    && visMass_muTau_1 > 80 && visMass_muTau_1<100 ) continue;
						
						
						if (taubyVLooseIsolationMVArun2017v2DBoldDMwLT2017->at(tauCand[0])==1  )
						  {
						    nVLooseTau+=event_weight;
						    h_tauPt_inc_vl->Fill(tauPt->at(tauCand[0]), event_weight);
						    h_tauM_vl->Fill(myTau.M());
						    if (tauDecayMode->at(tauCand[0])==0)
						      h_tauPt_dm0_vl->Fill(tauPt->at(tauCand[0]), event_weight);
						    if (tauDecayMode->at(tauCand[0])==1)
						      h_tauPt_dm1_vl->Fill(tauPt->at(tauCand[0]), event_weight);
						    if (tauDecayMode->at(tauCand[0])==10)
						      h_tauPt_dm10_vl->Fill(tauPt->at(tauCand[0]), event_weight);
						    if (tauDecayMode->at(tauCand[0])==11)
                                                      h_tauPt_dm11_vl->Fill(tauPt->at(tauCand[0]), event_weight);

						  }
						
						if (taubyTightIsolationMVArun2017v2DBoldDMwLT2017->at(tauCand[0])==1  )
						  {
						    nTightTau+=event_weight;
						    h_tauPt_inc_t->Fill(tauPt->at(tauCand[0]), event_weight);
						    h_tauM_t->Fill(myTau.M());
						    if (tauDecayMode->at(tauCand[0])==0)
						      h_tauPt_dm0_t->Fill(tauPt->at(tauCand[0]), event_weight);
						    if (tauDecayMode->at(tauCand[0])==1)
						      h_tauPt_dm1_t->Fill(tauPt->at(tauCand[0]), event_weight);
						    if (tauDecayMode->at(tauCand[0])==10)
						      h_tauPt_dm10_t->Fill(tauPt->at(tauCand[0]), event_weight);
						    if (tauDecayMode->at(tauCand[0])==11)
                                                      h_tauPt_dm11_t->Fill(tauPt->at(tauCand[0]), event_weight);

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
	
	
	
	
	  
	  //qcd sequential cuts
	  
	  //h_Events_level= new TH1F("Events_level","Events at each level of selection",20,0,20);
	  h_Events_level->SetBinContent(1, nInspected_genWeighted);
	  h_Events_level->SetBinContent(2, nMETFiltersPassed);
          h_Events_level->SetBinContent(3, nSingleTrgPassed);          
	  h_Events_level->SetBinContent(4, nGoodMuonPassed);
	  h_Events_level->SetBinContent(5, nGoodMuonPassed_2);	  
          h_Events_level->SetBinContent(6, nGoodMuMuPassed);
          h_Events_level->SetBinContent(7, nDeltaRPassed);
	  h_Events_level->SetBinContent(8, nPassedBjetVeto);
	  h_Events_level->SetBinContent(9, nPassedVisibleMasscut); 
	  h_Events_level->SetBinContent(10, nGoodTauPassed); 
	  // 
          h_Cutflow->SetBinContent(1, nGoodMuonPassed);
          h_Cutflow->SetBinContent(2, nGoodMuonPassed_2);
          h_Cutflow->SetBinContent(3, nGoodMuMuPassed);
          h_Cutflow->SetBinContent(4, nDeltaRPassed);
	  h_Cutflow->SetBinContent(6, nPassedBjetVeto);
	
	  tree->Fill();
	  
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
	
	// //g_tauFF_dm0->TGraphErrors(h_tauFF_dm0);
	// TGraphAsymmErrors* g_tauFF_dm0 = new TGraphAsymmErrors();
	// g_tauFF_dm0->Divide(h_tauPt_dm0_t,h_tauPt_dm0_vl);
	// g_tauFF_dm0->GetXaxis()->SetTitle("X title");
	// g_tauFF_dm0->GetYaxis()->SetTitle("Y title");
	// g_tauFF_dm0->SetTitle("tau fake factor dm0");
	// g_tauFF_dm0->Write("g_tauFF_dm0");


        // TGraphAsymmErrors* g_tauFF_dm1 = new TGraphAsymmErrors();
        // g_tauFF_dm1->Divide(h_tauPt_dm1_t,h_tauPt_dm1_vl);
        // g_tauFF_dm1->GetXaxis()->SetTitle("X title");
        // g_tauFF_dm1->GetYaxis()->SetTitle("Y title");
        // g_tauFF_dm1->SetTitle("tau fake factor dm1");
        // g_tauFF_dm1->Write("g_tauFF_dm1");

        // TGraphAsymmErrors* g_tauFF_dm10 = new TGraphAsymmErrors();
        // g_tauFF_dm10->Divide(h_tauPt_dm10_t,h_tauPt_dm10_vl);
        // g_tauFF_dm10->GetXaxis()->SetTitle("X title");
        // g_tauFF_dm10->GetYaxis()->SetTitle("Y title");
        // g_tauFF_dm10->SetTitle("tau fake factor dm10");
        // g_tauFF_dm10->Write("g_tauFF_dm10");





	std::cout.setf( std::ios::fixed, std:: ios::floatfield );
	if((nentriesToCheck-1)%reportEvery != 0)
		std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
	sw.Stop();
	std::cout<<"All events checked."<<std::endl;
	std::cout<<"*******************************************"<<std::endl;
	std::cout<<"******************Jithin's original*************************"<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"Initial entries "<<numberOfEvents<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"Inspected genWeightd "<<nInspected_genWeighted<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"METFiltersPassed "<<nMETFiltersPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"SingleTrgPassed "<<nSingleTrgPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"GoodMuonPassed "<<nGoodMuonPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"Good Muon 2 Passed "<<nGoodMuonPassed_2<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"opp charge "<<nGoodMuMuPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"DeltaRPassed "<<nDeltaRPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"PassedBjetVeto "<<nPassedBjetVeto<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"PassedVisibleMasscut "<<nPassedVisibleMasscut<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"nGoodTauPassed " << nGoodTauPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"nTightTauPassed " << nTightTau<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"nVLooseTauPassed " << nVLooseTau<<std::endl;

	std::cout<<"*******************************************"<<std::endl;
	std::cout<<"*******************************************"<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"Number of events inspected: " << nInspected <<std::endl;
	std::cout<<std::setw(20) <<std::right << "Number of events inspected (minus negative gen. weights): " << nInspected_genWeighted << std::endl; 
	std::cout<<std::setw(20) <<std::right << " (net weight with SF): " << netWeight << std::endl;
	
	std::cout<<"*******************************************"<<std::endl;
	std::cout<<std::setprecision(3)<<"  Real Time = " << sw.RealTime() << ",    Cpu Time = "<<sw.CpuTime()<< std::endl;
}

void Analyzer_mutau_data::BookHistos(const char* file2)
{
	fileName = new TFile(file2, "RECREATE");
	tree = new TTree("ADD","ADD");
	tree->Branch("event_","std::vector<unsigned int>",&event_);
	tree->Branch("event_info","std::vector<double>",&event_info);
	fileName->cd();

	Float_t PtBins[21]={0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,110, 120,130, 140,150, 160,170, 180,190, 200.};	
	Float_t MetBins[15]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180., 200, 300., 400., 600.0,800.0};
	Float_t TrMassBins[16]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180., 200, 300., 400., 600.0,800.0, 1000.0};

	/*	//Set up the histos to be filled with method fillHistos
	for(int i=0; i<10; i++)
	  {
		char ptbins[100];
		sprintf(ptbins, "_%d", i);
		std::string histname(ptbins);
		h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(), "nVtx",40,0,40);h_nVtx[i]->Sumw2();
		h_nEvents[i] = new TH1F(("nEvents"+histname).c_str(), "nEvents",3,0,3);h_nEvents[i]->Sumw2();
                //h_genWeight[i] = new TH1F(("genWeight"+histname).c_str(), "genWeight",10, -10.0, 10.0);h_genWeight[i]->Sumw2();
                //h_genHT[i] = new TH1F(("genHT"+histname).c_str(), "genHT",20,PtBins);h_genHT[i]->Sumw2();
		//**************  electrons  **************
		h_muon_En[i] = new TH1F(("Muon_En"+histname).c_str(), "Muon_En",20,PtBins);h_muon_En[i]->Sumw2();
		h_muon_Pt[i] = new TH1F(("Muon_Pt"+histname).c_str(), "Muon_Pt",20,PtBins);h_muon_Pt[i]->Sumw2();
		h_muon_eta[i] = new TH1F(("Muon_eta"+histname).c_str(), "Muon_eta",20,-3.0, 3.0);h_muon_eta[i]->Sumw2();
		h_muon_SCEta[i] = new TH1F(("Muon_SCeta"+histname).c_str(), "Muon_SCeta",20, -3.0, 3.0);h_muon_SCEta[i]->Sumw2();
		h_muon_phi[i] = new TH1F(("Muon_phi"+histname).c_str(), "Muon_phi", 21,-3.14,3.14);h_muon_phi[i]->Sumw2();
		h_muon_SCPhi[i] = new TH1F(("muon_SCphi"+histname).c_str(), "muon_SCphi", 21,-3.14,3.14);h_muon_SCPhi[i]->Sumw2();
		h_muon_IDbit[i] = new TH1F(("muon_ID_bit"+histname).c_str(), "muon_ID_bit",8,0,8);h_muon_IDbit[i]->Sumw2();
		//**************  tau  **************
		h_tau_En[i] = new TH1F(("Tau_En"+histname).c_str(), "Tau_En",20,PtBins);h_tau_En[i]->Sumw2();
		h_tau_Pt[i] = new TH1F(("Tau_Pt"+histname).c_str(), "Tau_Pt",20,PtBins);h_tau_Pt[i]->Sumw2();
		h_tau_eta[i] = new TH1F(("Tau_eta"+histname).c_str(), "Tau_eta",20,-3.0, 3.0);h_tau_eta[i]->Sumw2();

		h_tau_phi[i] = new TH1F(("Tau_phi"+histname).c_str(), "Tau_phi", 21,-3.14,3.14);h_tau_phi[i]->Sumw2();
		//**************  Met  **************
		h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "pfMET",14,MetBins);h_pfMET[i]->Sumw2();

		h_dPhi[i] = new TH1F(("h_dPhi"+histname).c_str(),"h_dPhi",20,0,3.15);h_dPhi[i]->Sumw2();
                h_dR[i] = new TH1F(("h_dR"+histname).c_str(),"h_dR",20,0,3.15);h_dR[i]->Sumw2();

		h_nJet[i] = new TH1F(("nJet"+histname).c_str(), "nJet",20,0,20);h_nJet[i]->Sumw2();
		h_leadingJetPt[i] = new TH1F(("leadingJetPt"+histname).c_str(),"leadingJetPt",30,20,1000);h_leadingJetPt[i]->Sumw2();

		h_leadingJetEta[i] = new TH1F(("h_leadingJetEta"+histname).c_str(),"h_leadingJetEta",40,-1.4442,1.4442);h_leadingJetEta[i]->Sumw2();

		//**************  rest  **************

		h_Mt[i]= new TH1F(("Mt"+histname).c_str(),"MT",15, TrMassBins);h_Mt[i]->Sumw2();
		h_VisibleMass[i]= new TH1F(("VisibleMass"+histname).c_str(),"VisibleMass",40,0,200);h_VisibleMass[i]->Sumw2();
		h_totTrMass[i]= new TH1F(("TotalTrMass"+histname).c_str(),"TotalTrMass",15,TrMassBins);h_totTrMass[i]->Sumw2();
		h_HiggsPt[i]= new TH1F(("HiggsPt"+histname).c_str(),"HiggsPt",20,PtBins);h_HiggsPt[i]->Sumw2();
		h_muonIso[i]= new TH1F(("Muon_iso"+histname).c_str(),"Muon_iso",40,0.0,1.0);h_muonIso[i]->Sumw2();
		h_tauIso[i]= new TH1F(("Tau iso"+histname).c_str(),"Tau_iso",12,0.0,1.2);h_tauIso[i]->Sumw2();
		h_tauMass[i]= new TH1F(("Tau_mass"+histname).c_str(),"Tau_mass", 18, 0.0, 1.8);h_tauMass[i]->Sumw2();
		h_tauDecayMode[i]= new TH1F(("Tau_Decay_Mode"+histname).c_str(),"Tau_Decay_Mode", 11, 0.0, 11.0);h_tauDecayMode[i]->Sumw2();


	}
	*/
}

//Fill the sequential histos at a particular spot in the sequence
void Analyzer_mutau_data::fillHistos(int histoNumber, double event_weight,int muIndex, int tauIndex)
{
  /*	//*********** fill electrons  ***********
		h_muon_En[histoNumber]->Fill((muEn->at(muIndex)),event_weight);
		h_muon_Pt[histoNumber]->Fill((muPt->at(muIndex)),event_weight);		
		h_muon_eta[histoNumber]->Fill(muEta->at(muIndex),event_weight);
		//h_muon_SCEta[histoNumber]->Fill(eleSCEta->at(muIndex),event_weight);
		h_muon_phi[histoNumber]->Fill(muPhi->at(muIndex),event_weight);
		//h_muon_SCPhi[histoNumber]->Fill(eleSCPhi->at(muIndex),event_weight);
		h_muon_IDbit[histoNumber]->Fill( muIDbit->at(muIndex)>>3&1,event_weight);
	//*********** fill taus  ***********

		h_tau_En[histoNumber]->Fill((tauEnergy->at(tauIndex)),event_weight);
		h_tau_Pt[histoNumber]->Fill((tauPt->at(tauIndex)),event_weight);
		h_tau_eta[histoNumber]->Fill(tauEta->at(tauIndex),event_weight);
		h_tau_phi[histoNumber]->Fill(tauPhi->at(tauIndex),event_weight);
		
	//*********** fill met  ***********		
		h_pfMET[histoNumber]->Fill(pfMET,event_weight);
		double dPhi_mutau = DeltaPhi(muPhi->at(muIndex),tauPhi->at(tauIndex));
		h_dPhi[histoNumber]->Fill(dPhi_mutau,event_weight);
		double dR_mutau = dR(muIndex,tauIndex);
		h_dR[histoNumber]->Fill(dR_mutau,event_weight);


	//*********** fill rest  ***********

		float mT_muMet = TMass_F((muPt->at(muIndex)),(muPhi->at(muIndex)),pfMET,pfMETPhi  );
		h_Mt[histoNumber]->Fill(mT_muMet,event_weight);
		TLorentzVector myTau; 
		myTau.SetPtEtaPhiE(tauPt->at(tauIndex),tauEta->at(tauIndex),tauPhi->at(tauIndex), tauEnergy->at(tauIndex));		
		TLorentzVector myMu; 
		myMu.SetPtEtaPhiE(muPt->at(muIndex),muEta->at(muIndex),muPhi->at(muIndex), muEn->at(muIndex));
		double visMass_mutau = VisMass_F(myTau, myMu);
		h_VisibleMass[histoNumber]->Fill(visMass_mutau,event_weight);
		TLorentzVector myMet;
                myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
		float tot_tr_mass = TotTMass_F(myMu, myTau, myMet );
		h_totTrMass[histoNumber]->Fill(tot_tr_mass,event_weight);


		double HiggsPt = pTvecsum_F(muPt->at(muIndex),tauPt->at(tauIndex),muPhi->at(muIndex),tauPhi->at(tauIndex) );
		h_HiggsPt[histoNumber]->Fill(HiggsPt,event_weight);
		h_nVtx[histoNumber]->Fill(nVtx,event_weight);
		h_nEvents[histoNumber]->Fill(isData,event_weight);
		
		h_nJet[histoNumber]->Fill(nJet ,event_weight);
		//		if(jetsize>0){
		//h_leadingJetPt[histoNumber]->Fill(jetPt->at(),event_weight);
		//h_leadingJetPt_300[histoNumber]->Fill(jetPt->at(),event_weight);
		//h_leadingJetEta[histoNumber]->Fill(jetEta->at(),event_weight);
		//}
		//h_genHT[histoNumber]->Fill(genHT,event_weight);
		//h_genWeight[histoNumber]->Fill(genWeight,event_weight);
		float rel_muIso = ( muPFChIso->at(muIndex) + max( muPFNeuIso->at(muIndex) + muPFPhoIso->at(muIndex) - 0.5 *muPFPUIso->at(muIndex) , 0.0 )) / (muPt->at(muIndex));
		h_muonIso[histoNumber]->Fill(rel_muIso,event_weight);
		h_tauIso[histoNumber]->Fill(tauByTightIsolationMVArun2v1DBoldDMwLT->at(tauIndex),event_weight);
		h_tauMass[histoNumber]->Fill(tauMass->at(tauIndex),event_weight);
		h_tauDecayMode[histoNumber]->Fill(tauDecayMode->at(tauIndex),event_weight);

  */
}




//---------------------------------------------------                                                                                                                                
// get a electron candiate based on pt eta and isolation                                                                                                                               
//----------------------------------------------------                                                                                                                               

std::vector<int> Analyzer_mutau_data::getMuCand_1(double muPtCut, double muEtaCut){
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
     
      float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));

      if( relMuIso < 0.15 ) relative_iso = true;
      if(muonId && kinematic && relative_iso){
	tmpCand.push_back(iMu);
      }                                                                                      
    }                                                                                       
  return tmpCand;
}
std::vector<int> Analyzer_mutau_data::getMuCand_2(double muPtCut, double muEtaCut, int skip_i){
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

      float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));

      if( relMuIso < 0.15 ) relative_iso = true;
      if(muonId && kinematic && relative_iso){
        tmpCand.push_back(iMu);
      }
    }
  return tmpCand;
}


std::vector<int> Analyzer_mutau_data::getTauCand(double tauPtCut, double tauEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over taus      
  for(int iTau=0;iTau<nTau;iTau++)
    {
      bool kinematic = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;

      if( tauPt->at(iTau) > tauPtCut  && fabs( tauEta->at(iTau))< tauEtaCut && fabs(taudz->at(iTau))<0.2 )kinematic = true;
      if( tauByTightMuonRejection3->at(iTau) == 1 && tauByMVA6LooseElectronRejection->at(iTau)==1) tauId = true;
      //if( (taubyTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau)==1)) tauIsolation = true;
      //if( tauDecayMode->at(iTau)==0 || tauDecayMode->at(iTau)==1 || tauDecayMode->at(iTau)==10 ) decayModeCut = true;
      
      if(kinematic == true && tauId==true )
	//if(tauId==true && kinematic==true && decayModeCut==true)
	{
	  tmpCand.push_back(iTau);
	}                                                           
    }                                                                                       
  return tmpCand;
}


std::vector<int> Analyzer_mutau_data::getASRTauCand(double tauPtCut, double tauEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over taus
  for(int iTau=0;iTau<nTau;iTau++)
    {
      bool kinematic = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;

      if( tauPt->at(iTau) > tauPtCut  && fabs( tauEta->at(iTau))< tauEtaCut && taudz->at(iTau)<0.2 )kinematic = true;
      if( tauByTightMuonRejection3->at(iTau) == 1 && tauByMVA6VLooseElectronRejection->at(iTau) == 1) tauId = true;
      if( (taubyTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau)!=1 && taubyVLooseIsolationMVArun2017v2DBoldDMwLT2017->at(iTau)==1)) tauIsolation = true;
      if( tauDecayMode->at(iTau)==0 || tauDecayMode->at(iTau)==1 || tauDecayMode->at(iTau)==10) decayModeCut = true;

      if(tauId==true && kinematic==true && tauIsolation==true && decayModeCut==true)
        //if(tauId==true && kinematic==true && decayModeCut==true)
        {
          tmpCand.push_back(iTau);
        }
    }
  return tmpCand;
}


bool Analyzer_mutau_data::thirdLeptonVeto(int muon_1, int muon_2)
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
       if( muIDbit->at(iMu)>>15&1==1) muonId =true;
       bool relative_iso = false;

       float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
       if( relMuIso < 0.3 ) relative_iso = true;


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

double Analyzer_mutau_data::dR(int mu_1_index, int mu_2_index )
{
  double deltaeta = abs(muEta->at(mu_1_index) - muEta->at(mu_2_index));
  double muon_1_Phi = muPhi->at(mu_1_index);
  double muon_2_Phi = muPhi->at(mu_2_index);

  double deltaphi = DeltaPhi(muon_1_Phi, muon_2_Phi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}



double Analyzer_mutau_data::dR_muTau(int mu_2_index, int tau_index )
{
  double deltaeta = abs(muEta->at(mu_2_index) - tauEta->at(tau_index));
  double muon_2_Phi = muPhi->at(mu_2_index);
  double tau_Phi = tauPhi->at(tau_index);

  double deltaphi = DeltaPhi(muon_2_Phi, tau_Phi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

double Analyzer_mutau_data::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
}

float Analyzer_mutau_data::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
    return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

float Analyzer_mutau_data::TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float totalTMass = (a + b+ met).M();
  return totalTMass;
}


float Analyzer_mutau_data::VisMass_F(TLorentzVector a, TLorentzVector b){
  float visibleMass = (a + b).M();
  return visibleMass;
}

float Analyzer_mutau_data::pTvecsum_F(float pt1, float pt2, float phi1, float phi2) {
  float pt_vecSum = sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2), 2) + pow(pt1*sin(phi1) + pt2*sin(phi2), 2));
  return pt_vecSum;
}

bool Analyzer_mutau_data::passBjetVeto()
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
