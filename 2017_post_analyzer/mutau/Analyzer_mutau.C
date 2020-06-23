////Analyzer_mutau.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom Analyzer_mutau analyze
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
#define Analyzer_mutau_cxx
#include "Analyzer_mutau.h"
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
#include <set>
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
#include <TLorentzVector.h>

using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{

  std::string SampleName = argv[5];

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
  Analyzer_mutau t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery, SampleName);
  return 0;
}

void Analyzer_mutau::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
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
   double nHLTPassed,n_eventWeight, nSingleTrgPassed, nGoodMuonPassed, nElectronPtPassed, nGoodTauPassed, nTauPtPassed,numberOfEvents,nMETPassed, nDPhiPassed, nqcdden,nqcdnum, nMETFiltersPassed,nLeptonVetoPassed,nPassedBjetVeto,nNoisyCrystals,nDPhiJetMETPassed, nGoodMuTauPassed,nDeltaRPassed,nPassedThirdLepVeto,nMtPassed,  nPassedHiggsPtcut, nPassedVisibleMasscut, nPassedMETcut, nFinal_afterSelections, nGoodMuonPassed_qcd, nGoodTauPassed_qcd, nGoodMuTauPassed_qcd, nDeltaRPassed_qcd, nPassedThirdLepVeto_qcd,nPassedBjetVeto_qcd, nPassedHiggsPtcut_qcd, nPassedVisibleMasscut_qcd, nPassedMETcut_qcd, nFinal_afterSelections_qcd ;
   nHLTPassed = n_eventWeight = nSingleTrgPassed = nGoodMuonPassed = nElectronPtPassed = nGoodTauPassed = nTauPtPassed= numberOfEvents = nMETPassed = nDPhiPassed = nqcdden= nqcdnum=nMETFiltersPassed= nLeptonVetoPassed=nPassedBjetVeto=nNoisyCrystals=nDPhiJetMETPassed= nGoodMuTauPassed = nDeltaRPassed= nPassedThirdLepVeto=nMtPassed=nPassedHiggsPtcut=nPassedVisibleMasscut=nPassedMETcut=nFinal_afterSelections=nGoodMuonPassed_qcd=nGoodTauPassed_qcd=nGoodMuTauPassed_qcd=nDeltaRPassed_qcd=nPassedThirdLepVeto_qcd=nPassedBjetVeto_qcd=nPassedHiggsPtcut_qcd=nPassedVisibleMasscut_qcd=nPassedMETcut_qcd=nFinal_afterSelections_qcd=0;

   std::vector<int> muCand;
   muCand.clear();
   std::vector<int> tauCand;
   tauCand.clear();


   TString sample = TString(SampleName);

   Float_t met_bins[16]={0.0, 15, 30, 45, 60, 90, 105, 120, 135, 150,175., 190., 250., 400., 600.0,800.0};
   TH1F* h_Events_level= new TH1F("Events_level","Events at each level of selection",15,0,30);
   TH1F* h_Cutflow= new TH1F("Cutflow","Events at each level of selection",6,0,12);
   TH1F* h_metfilter= new TH1F("metfilter","metfilter values",1000,0,1000); 
   TH1F* h_Events_level_qcd= new TH1F("Events_level_qcd","Events at each level of selection",15,0,30);
   TH1F* h_Cutflow_qcd= new TH1F("Cutflow_qcd","Events at each level of selection",6,0,12);

       
   TH1F* h_MET_0= new TH1F("MET_0","MET before any selection",15,met_bins );  h_MET_0->Sumw2();
   TH1F* h_MET_1= new TH1F("MET_1","MET after met filters",15,met_bins ); h_MET_1->Sumw2();
   TH1F* h_MET_2= new TH1F("MET_2","MET after single ele trg",15,met_bins ); h_MET_2->Sumw2();
   TH1F* h_MET_3= new TH1F("MET_3","MET after ele selection",15,met_bins ); h_MET_3->Sumw2();
   TH1F* h_MET_4= new TH1F("MET_4","MET after tau selection",15,met_bins ); h_MET_4->Sumw2();


   //   TH1F* h_mupt_0= new TH1F("mupt_0","mupt_0",40,0,400 ); h_mupt_0->Sumw2();
   //TH1F* h_mupt_1= new TH1F("mupt_1","mupt_1",40,0,400 ); h_mupt_1->Sumw2();
   //TH1F* h_mupt_2= new TH1F("mupt_2","mupt_2",40,0,400 ); h_mupt_2->Sumw2();


   //TH1F* h_genWeight = new TH1F((sample+"genWeight").c_str(), "genWeight",3000, -30000.0, 30000.0);h_genWeight[i]->Sumw2();
   //TH1F* h_genHT = new TH1F((sample+"genHT").c_str(), "genHT",15,PtBins);h_genHT[i]->Sumw2();
   TH1F* h_isData= new TH1F("isData","isData",4,0.0,2.0 );  h_isData->Sumw2();
   TH1F* h_insEvents= new TH1F("insEvents","insEvents",2,0.0,2.0 );  h_insEvents->Sumw2();

   int iphi = 41;
   int ieta = 5;


   TFile *f_muIDSF=new TFile("RunBCDEF_SF_ID.root");
   TH2F *h_muIDSF=(TH2F*) f_muIDSF->Get("NUM_TightID_DEN_genTracks_pt_abseta");

   TFile *f_muIsoSF=new TFile("RunBCDEF_SF_ISO.root");
   TH2F *h_muIsoSF=(TH2F*) f_muIsoSF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

   TFile *f_muTrgSF=new TFile("EfficienciesAndSF_RunBtoF_Nov17Nov2017.root");
   TH2F *h_muTrgSF=(TH2F*) f_muTrgSF->Get("IsoMu27_PtEtaBins/pt_abseta_ratio");


   //   TFile *f_kfactors = new TFile("kfactors.root");
   //TH1D *ewkCorrection = (TH1D*)f_kfactors->Get("EWKcorr/W");
   //TH1D *NNLOCorrection = (TH1D*)f_kfactors->Get("WJets_LO/inv_pt");

   //std::cout<<"Pont 1" << endl;
   //bool debug=true;
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
	  double inspected_event_weight = 1.0; 
	  fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight/fabs(genWeight) : inspected_event_weight = 0.0;
	  nInspected_genWeighted += inspected_event_weight;  
	  nInspected += 1; 
	  h_insEvents->SetBinContent(1, nInspected_genWeighted);
	  //=1.0 for real data
	  double event_weight=1.0;
	  double sf_tauID = 1.0; 
	  double sf_IsoEff = 1.0; 
	  double sf_muTrg = 1.0;
	  double sf_muID = 1.0;

	  /*
	  double EWK_corrected_weight=1.0;
	  double NNLO_weight = 1.0;
	  double kfactor = 1.0;

	  int bosonPID;
	  double bosonPt;
	  bool Wfound = false;
	  //check which mc particle is W boson                                                                                               
	  for(int i=0; i<nMC;i++){
	    if((*mcPID)[i] == 24){
	      Wfound=true;
	      bosonPID = (*mcPID)[i];
	      bosonPt = (*mcPt)[i];
	    }
	  }
	  */
	  //std::cout<<"This works! P 1" <<endl;
	  muCand = getMuCand(30,2.4); // Electron pt, eta cuts
	  tauCand = getTauCand(20,2.3); // Tau pt, eta cuts
	  //particleCand = getParticleCand();
	  //	  if (run > 274421 && (event % 4) != 0) continue;
	  //fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;   
	  numberOfEvents+=event_weight;
	  
	  
	  
	  h_MET_0->Fill(pfMET, event_weight);
	  h_isData->Fill(isData, event_weight);
	  if(metFilters==0) 
	    {	    
	      fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	      nMETFiltersPassed+=event_weight;
	      h_MET_1->Fill(pfMET, event_weight);
	      if (debug==true)if (jentry%reportEvery == 0){ std::cout<<"This works! Met filters" <<endl;}

	      if(HLTEleMuX>>19&1 == 1 ) // Single muon trigger
		{
		  nSingleTrgPassed+=event_weight;
		  if (debug==true)if (jentry%1000 == 0) std::cout<<"This works! Single electron trigger" <<endl;
		  h_MET_2->Fill(pfMET, event_weight);


		  //###################### SIGNAL REGION ###############
		  if(muCand.size() >0) 
		    {
		      nGoodMuonPassed+=event_weight;
		      //std::cout<<"This works! muCand" <<endl;
		      //h_mupt_0->Fill(muPt->at(muCand[0]), event_weight); 
		      h_MET_3->Fill(pfMET, event_weight);
		      // if (jentry%1000 == 0) std::cout<<"This works! P 2" <<endl;
		      if( tauCand.size() >0 )
			{
			  //std::cout<<"This works! tauCand" <<endl;
			  nGoodTauPassed+=event_weight;
			  h_MET_4->Fill(pfMET, event_weight);
			  //h_mupt_1->Fill(muPt->at(muCand[0]), event_weight);
			  if (debug==true)std::cout<<"event_weight (PE) =  "<< event_weight<<std::endl;
			  
			  if (muPt->at(muCand[0]) <=120)
			    {
			      sf_muID = h_muIDSF->GetBinContent(h_muIDSF->GetXaxis()->FindBin(muPt->at(muCand[0])),h_muIDSF->GetYaxis()->FindBin(abs(muEta->at(muCand[0]))));
			      if (debug==true)std::cout<<"sf_muID =  "<< sf_muID<<std::endl;
			      sf_IsoEff = h_muIsoSF->GetBinContent(h_muIDSF->GetXaxis()->FindBin(muPt->at(muCand[0])),h_muIDSF->GetYaxis()->FindBin(abs(muEta->at(muCand[0]))));
			      if (debug==true)std::cout<<"sf_IsoEff =  "<< sf_IsoEff<<std::endl;
			    }
			  sf_muTrg = h_muTrgSF->GetBinContent(h_muIDSF->GetXaxis()->FindBin(muPt->at(muCand[0])),h_muIDSF->GetYaxis()->FindBin(abs(muEta->at(muCand[0]))));
			  if (debug==true)std::cout<<"sf_muTrg =  "<< sf_muTrg<<std::endl;
			  
			  event_weight=event_weight*sf_tauID*sf_muID*sf_IsoEff*sf_muTrg;
			  if (debug==true)std::cout<<"event_weight (AFTER SF) =  "<< event_weight<<std::endl;
			  
			  netWeight = event_weight;
			  //h_mupt_2->Fill(muPt->at(muCand[0]), event_weight);
			  if( (muCharge->at(muCand[0]))*(tauCharge->at(tauCand[0])) < 0 ) // opposite cgharge condition
			    {
			      if (debug==true)std::cout<<"This works! Charge criteria" <<endl;
			      //h_mupt_2->Fill(muPt->at(muCand[0]), event_weight);
			      nGoodMuTauPassed+=event_weight;
			      fillHistos(0,event_weight,muCand[0], tauCand[0]);
			      double deltaR = dR(muCand[0],tauCand[0]);
			      if(deltaR > 0.3 )
				{
				  if (debug==true)std::cout<<"This works! dR cut" <<endl;
				  nDeltaRPassed+=event_weight;
				  fillHistos(1,event_weight,muCand[0], tauCand[0]);
				  
				  if( thirdLeptonVeto()==true )
				    {
				      nPassedThirdLepVeto+=event_weight;
				      if (debug==true)std::cout<<"This works! third lepon veto" <<endl;
				      fillHistos(2,event_weight,muCand[0], tauCand[0]);
				      if( passBjetVeto() == true)
					{
					  nPassedBjetVeto+=event_weight;
					  if (debug==true)if (jentry%5 == 0) std::cout<<"This works! bJetVeto" <<endl;
					  fillHistos(3,event_weight,muCand[0], tauCand[0]);
					  if (debug==true)if (jentry%5 == 0) std::cout<<"This works! bJetVeto **** after fill" <<endl;
					  
					  
					  float muMet_trMass = TMass_F((muPt->at(muCand[0])),(muPhi->at(muCand[0])),pfMET,pfMETPhi  );
					  //if( (muMet_trMass>50) )
					    {
					      fillHistos(4,event_weight,muCand[0], tauCand[0]);
					      nMtPassed+=event_weight;
					      double HiggsPt = pTvecsum_F(muPt->at(muCand[0]),tauPt->at(tauCand[0]),muPhi->at(muCand[0]),tauPhi->at(tauCand[0]));
					      if( (HiggsPt>65) ) 
						{
						  nPassedHiggsPtcut+=event_weight;
						  fillHistos(5,event_weight,muCand[0], tauCand[0]);
						  
						  TLorentzVector myTau;
						  myTau.SetPtEtaPhiE(tauPt->at(tauCand[0]),tauEta->at(tauCand[0]),tauPhi->at(tauCand[0]), tauEnergy->at(tauCand[0]));
						  TLorentzVector myMu;
						  myMu.SetPtEtaPhiE(muPt->at(muCand[0]),muEta->at(muCand[0]),muPhi->at(muCand[0]), muEn->at(muCand[0]));
						  double visMass_etau = VisMass_F(myTau, myMu);
						  
						  if((visMass_etau < 125) ) 
						    {
						      nPassedVisibleMasscut+=event_weight;
						      fillHistos(6,event_weight,muCand[0], tauCand[0]);
						      
						      
						      if( (pfMET > 105) ) 
							{
							  nPassedMETcut+=event_weight;
							  fillHistos(7,event_weight,muCand[0], tauCand[0]);
							  
							  if (debug==true)if (jentry%5 == 0) std::cout<<"******************* This works! All selections passed ***************" <<endl;
							  nFinal_afterSelections+=event_weight;
							  
							  TLorentzVector myTau_1; 
							  myTau_1.SetPtEtaPhiE(tauPt->at(tauCand[0]),tauEta->at(tauCand[0]),tauPhi->at(tauCand[0]), tauEnergy->at(tauCand[0]));		
							  TLorentzVector myMu_1; 
							  myMu_1.SetPtEtaPhiE(muPt->at(muCand[0]),muEta->at(muCand[0]),muPhi->at(muCand[0]), muEn->at(muCand[0]));
							  TLorentzVector myMet_1;
							  myMet_1.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
							  float tot_tr_mass = TotTMass_F(myMu_1, myTau_1, myMet_1 );
							  
							  if (tot_tr_mass > 260)
							    {
							      fillHistos(8,event_weight,muCand[0], tauCand[0]);
							      
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
	    
	
	  //qcd sequential cuts
	  
	  //h_Events_level= new TH1F("Events_level","Events at each level of selection",20,0,20);
	  h_Events_level->SetBinContent(1, nInspected_genWeighted);
	  h_Events_level->SetBinContent(2, nMETFiltersPassed);
          h_Events_level->SetBinContent(3, nSingleTrgPassed);          
	  h_Events_level->SetBinContent(4, nGoodMuonPassed);
	  h_Events_level->SetBinContent(5, nGoodTauPassed);	  
          h_Events_level->SetBinContent(6, nGoodMuTauPassed);
          h_Events_level->SetBinContent(7, nDeltaRPassed);
          h_Events_level->SetBinContent(8, nPassedThirdLepVeto); 
	  h_Events_level->SetBinContent(9, nPassedBjetVeto);
	  h_Events_level->SetBinContent(10, nMtPassed);
	  h_Events_level->SetBinContent(11, nPassedHiggsPtcut); 
	  h_Events_level->SetBinContent(12, nPassedVisibleMasscut); 
	  h_Events_level->SetBinContent(13, nPassedMETcut); 
	  // 
          h_Cutflow->SetBinContent(1, nGoodMuonPassed);
          h_Cutflow->SetBinContent(2, nGoodTauPassed);
          h_Cutflow->SetBinContent(3, nGoodMuTauPassed);
          h_Cutflow->SetBinContent(4, nDeltaRPassed);
          h_Cutflow->SetBinContent(5, nPassedThirdLepVeto);
          h_Cutflow->SetBinContent(6, nPassedBjetVeto);
	
	  tree->Fill();
	  
	  if (jentry%reportEvery == 0)
	    {
	      std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
	    }
	}
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
	std::cout<<std::setw(20) <<std::right <<"GoodElectronPassed "<<nGoodMuonPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"GoodTauPassed "<<nGoodTauPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"opp charge "<<nGoodMuTauPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"DeltaRPassed "<<nDeltaRPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"PassedThirdLepVeto "<<nPassedThirdLepVeto<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"PassedBjetVeto "<<nPassedBjetVeto<<std::endl;
	std::cout<<"*******************************************"<<std::endl;
	std::cout<<"*******************************************"<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"Number of events inspected: " << nInspected <<std::endl;
	std::cout<<std::setw(20) <<std::right << "Number of events inspected (minus negative gen. weights): " << nInspected_genWeighted << std::endl; 
	std::cout<<std::setw(20) <<std::right << " (net weight with SF): " << netWeight << std::endl;
	std::cout<<"*******************************************"<<std::endl;

}

void Analyzer_mutau::BookHistos(const char* file2)
{
	fileName = new TFile(file2, "RECREATE");
	tree = new TTree("ADD","ADD");
	tree->Branch("event_","std::vector<unsigned int>",&event_);
	tree->Branch("event_info","std::vector<double>",&event_info);
	fileName->cd();

	Float_t PtBins[21]={0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,110, 120,130, 140,150, 160,170, 180,190, 200.};	
	Float_t MetBins[15]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180., 200, 300., 400., 600.0,800.0};
	Float_t TrMassBins[16]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180., 200, 300., 400., 600.0,800.0, 1000.0};

	//Set up the histos to be filled with method fillHistos
	for(int i=0; i<10; i++)
	  {
		char ptbins[100];
		sprintf(ptbins, "_%d", i);
		std::string histname(ptbins);
		h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(), "nVtx",40,0,40);h_nVtx[i]->Sumw2();
		h_nEvents[i] = new TH1F(("nEvents"+histname).c_str(), "nEvents",3,0,3);h_nEvents[i]->Sumw2();
                h_genWeight[i] = new TH1F(("genWeight"+histname).c_str(), "genWeight",10, -10.0, 10.0);h_genWeight[i]->Sumw2();
                h_genHT[i] = new TH1F(("genHT"+histname).c_str(), "genHT",20,PtBins);h_genHT[i]->Sumw2();
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
}

//Fill the sequential histos at a particular spot in the sequence
void Analyzer_mutau::fillHistos(int histoNumber, double event_weight,int muIndex, int tauIndex)
{
	//*********** fill electrons  ***********
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
		h_genHT[histoNumber]->Fill(genHT,event_weight);
		h_genWeight[histoNumber]->Fill(genWeight,event_weight);
		float rel_muIso = ( muPFChIso->at(muIndex) + max( muPFNeuIso->at(muIndex) + muPFPhoIso->at(muIndex) - 0.5 *muPFPUIso->at(muIndex) , 0.0 )) / (muPt->at(muIndex));
		h_muonIso[histoNumber]->Fill(rel_muIso,event_weight);
		h_tauIso[histoNumber]->Fill(tauByTightIsolationMVArun2v1DBoldDMwLT->at(tauIndex),event_weight);
		h_tauMass[histoNumber]->Fill(tauMass->at(tauIndex),event_weight);
		h_tauDecayMode[histoNumber]->Fill(tauDecayMode->at(tauIndex),event_weight);

}




//---------------------------------------------------                                                                                                                                
// get a electron candiate based on pt eta and isolation                                                                                                                               
//----------------------------------------------------                                                                                                                               

std::vector<int> Analyzer_mutau::getMuCand(double muPtCut, double muEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over muons                                                                     
  for(int iMu=0;iMu<nMu;iMu++)
    {
      bool kinematic = false;
      if( (*muPt)[iMu] > muPtCut  && fabs((*muEta)[iMu])< muEtaCut && muDz->at(iMu) < 0.2 && muD0->at(iMu)<0.045) kinematic = true;
      bool muonId =false;
      if( muIDbit->at(iMu)>>3&1==1) muonId =true;
      bool relative_iso = false;
     
      float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));

      if( relMuIso < 0.15 ) relative_iso = true;
      if(muonId && kinematic && relative_iso){
	tmpCand.push_back(iMu);
      }                                                                                      
    }                                                                                       
  return tmpCand;
}


std::vector<int> Analyzer_mutau::getTauCand(double tauPtCut, double tauEtaCut){
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
      if( (taubyTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau)==1)) tauIsolation = true;
      if( tauDecayMode->at(iTau)==0 || tauDecayMode->at(iTau)==1 || tauDecayMode->at(iTau)==10 ) decayModeCut = true;
      
      if(tauId==true && kinematic==true && tauIsolation==true && decayModeCut==true)
	//if(tauId==true && kinematic==true && decayModeCut==true)
	{
	  tmpCand.push_back(iTau);
	}                                                           
    }                                                                                       
  return tmpCand;
}
std::vector<int> Analyzer_mutau::getASRTauCand(double tauPtCut, double tauEtaCut){
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


bool Analyzer_mutau::thirdLeptonVeto(){

  //  std::vector<int> tmpCand;
  //tmpCand.clear();
  int tmpCand = 0;
  bool thirdLep = false;
  bool thirdLepVeto=true;
   //Loop over muons                                                                                                                       

   for(int iEle=0; iEle < nEle;iEle++)
     {
       bool kinematic = false;
       if( (*elePt)[iEle] > 10.0  && fabs((*eleEta)[iEle])< 2.5 && (*eleD0)[iEle] < 0.045 && (*eleDz)[iEle] < 0.2 ) kinematic = true;
       bool electronId =false;
       if( eleIDbit->at(iEle)>>0&1==1) electronId =true;
       bool relative_iso = false;

       float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
       if( relEleIso < 0.3 ) relative_iso = true;
       if(electronId==true && kinematic==true && relative_iso==true){
	 //tmpCand.push_back(iEle);
	 tmpCand++;
       }                                                                                                                                                    
     }                                                                                                                      
   if(tmpCand > 0){ thirdLep = true; thirdLepVeto=false;}
   //   if( thirdLep = true ) thirdLepVeto=false;
   //if( tmpCand = 0 ) thirdLepVeto=true;

   return thirdLepVeto;

}




//Veto failed if a electron is found that passes Loose Muon ID, Loose Muon Isolation, and elePtcut, and does not overlap the candidate photon within dR of 0.5                                                                                                                                                                    

double Analyzer_mutau::dR(int mu_index, int tau_index)
{
  double deltaeta = abs(muEta->at(mu_index) - tauEta->at(tau_index));
  double muon_Phi = muPhi->at(mu_index);
  double tau_Phi = tauPhi->at(tau_index);

  double deltaphi = DeltaPhi(muon_Phi, tau_Phi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}



double Analyzer_mutau::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
}

float Analyzer_mutau::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
    return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

float Analyzer_mutau::TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float totalTMass = (a + b+ met).M();
  return totalTMass;
}


float Analyzer_mutau::VisMass_F(TLorentzVector a, TLorentzVector b){
  float visibleMass = (a + b).M();
  return visibleMass;
}

float Analyzer_mutau::pTvecsum_F(float pt1, float pt2, float phi1, float phi2) {
  float pt_vecSum = sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2), 2) + pow(pt1*sin(phi1) + pt2*sin(phi2), 2));
  return pt_vecSum;
}

bool Analyzer_mutau::passBjetVeto()
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
