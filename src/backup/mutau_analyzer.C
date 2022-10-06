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
//#include "roCorr_Run2_v3/RoccoR.cc"

#include "commonFunctions.h"
#include "fractions.C"

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
  
  mutau_analyzer t(argv[1],argv[2], isMC, SampleName);
  t.Loop(maxEvents,reportEvery, SampleName );
  //delete myMap1;
  cout<<" Outpt written to "<<outputfile<<endl;
  sw.Stop();
  sw.Print();
  return 0;
}

void mutau_analyzer::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
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
  if (fChain == 0) return;

  std::vector<int> muCand;        muCand.clear();
  std::vector<int> tauCand;       tauCand.clear();
  std::vector<int> aisrtauCand;   aisrtauCand.clear();
  TString sample = TString(SampleName);
  
  int nHiggs = 0;
  int nHToMuTau = 0;
  int found_mt = 0;
  int muCand_1=0; int muCand_2=0;int muCand_3=0;
  int tauCand_1=0; int tauCand_2=0;int tauCand_3=0;
  
  Double_t  Pt_Bins[26]={0.0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
  Double_t  Pt_Bins_highPt[21]={100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};

  //TH1F* h_cutflow=new TH1F("cutflow", "cutflow", 10, 0, 10); h_cutflow->Sumw2();
  TH1F* h_cutflow_n=new TH1F("cutflow_n", "cutflow_n", 8, 0, 8);h_cutflow_n->Sumw2();
  TH1F* h_cutflow_n_fr=new TH1F("cutflow_n_fr", "cutflow_n_fr", 8, 0, 8);h_cutflow_n_fr->Sumw2();
  TH1F* h_cutflow_n_dyll=new TH1F("cutflow_n_dyll", "cutflow_n_dyll", 8, 0, 8);h_cutflow_n_dyll->Sumw2();
  TH1F* h_cutflow_n_dyll_fr=new TH1F("cutflow_n_dyll_fr", "cutflow_n_dyll_fr", 8, 0, 8);h_cutflow_n_dyll_fr->Sumw2();
  //TH1F* h_cutflow_Htt=new TH1F("cutflow_Htt", "cutflow_Htt", 11, 0, 11); h_cutflow_Htt->Sumw2();

   
  // if(debug)cout<<" setting up other files ..."<<endl;
  //RoccoR  rc("sf_files/roCorr_Run2_v3/RoccoR2017.txt"); 
   
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
      if(debug) cout<<"event "<<jentry<<endl;
      muCand.clear(); 
      tauCand.clear();
     
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      double inspected_event_weight = 1.0; 
      if(is_MC)	 fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight/fabs(genWeight) : inspected_event_weight = 0.0;
      nInspected_genWeighted += inspected_event_weight;  
      nInspected += 1; 
      //h_insEvents->SetBinContent(1, nInspected_genWeighted);
      //=1.0 for real data
      double event_weight=1.0;
      double weight=1.0;
      double pileup_sf=1.0;
      double applySf=1.0;
      bool passSingleTriggerPaths=false;   bool passCrossTrigger=false;   int report_i=0;
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
       
      /////Trigger bit selection
      if(HLTEleMuX>>21&1 == 1 || HLTEleMuX>>60&1 == 1 )
	passSingleTriggerPaths=true;
      if( HLTTau>>0&1 == 1  )
	passCrossTrigger=true;
      ////
       
      ////// get muon, tau candidates from the list 
      muCand = getMuCand(20,2.1);
      tauCand = getTauCand(30,2.3);
      aisrtauCand = getAISRTauCand(30,2.3);


      if(!is_MC)
	event_weight=1.0;
      else	
	event_weight=weight;
      if(debug)cout<<"reco selections begin"<<endl;
       
      //muCand.clear(); tauCand.clear();
      //cout<<__LINE__<<endl;
      ////// DY Z-> ll signal region -  isolated begin
      //cout<<"entry: "<<jentry<<" , event : "<<event<<endl;
      if(metFilters==0 )
	{
	  if(debug)cout<<"metfilters selected"<<endl;
	  //cout<<__LINE__<<endl;
	  if(is_MC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	  nMETFiltersPassed_dyll+=event_weight;
	  if(debug)cout<<"genweight applied"<<endl;
	  if(   passSingleTriggerPaths || passCrossTrigger )
	    {
	      nSingleTrgPassed_dyll+=event_weight;
	      if(debug)cout<<"trigger selected"<<endl;
	      //muCand = getMuCand(20,2.4);  ///// muons selected 
	      if( muCand.size() >0 ) 
		{ 
		  nGoodMuonPassed_dyll+=event_weight;
		  if(debug)cout<<"this worked Line 443"<<endl;
		   
		  //tauCand = getTauCand(30,2.3);
		  if( tauCand.size()>0 ) 
		    {
		      //nGoodTauPassed_dyll+=event_weight;
		      if(debug)cout<<"this worked Line 305"<<endl;
		      setMyEleTau(muCand[0], tauCand[0]);
		      //printTabSeparated("in Z->ll", MuIndex, TauIndex, my_muP4.Pt(), my_tauP4.Pt(), my_metP4.Pt(), my_njets );
		      // if( passDiMuonVeto(MuIndex)==true 
		      // 	  && eVetoZTTp001dxyz(MuIndex, TauIndex)
		      // 	  && mVetoZTTp001dxyz(MuIndex, TauIndex)
		      // 	  ) Ztt_selector=true;
		      // else Ztt_selector=false;
		      
		      if ( TriggerSelection(my_muP4, my_tauP4) )
			{
			  //if(Ztt_selector) 
			  if (pass3rdLeptonVeto)
			    {
			      nGoodTauPassed_dyll+=event_weight;
			      if ( muCharge->at(MuIndex) * tau_Charge->at(TauIndex) < 0  
				   &&  (if_DY_Genmatching(MuIndex, TauIndex)==1 || if_DY_Genmatching(MuIndex, TauIndex)==2)  )
				{
				  nGoodMuTauPassed_dyll+=event_weight;
				  //makeTestPlot("e_dyll", 0,0,0,event_weight);
				   
				  if ( MatchTriggerFilter(MuIndex, TauIndex) )
				    {
				       
				      if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
				       
				      //my_muP4 = my_muP4*muRC_sf;
				      applySf=1.0;
				      if(is_MC)
					applySf=  getScaleFactors( my_muP4.Pt(),
								   my_tauP4.Pt(),
								   my_muP4.Eta(),
								   my_tauP4.Eta(),
								   tau_DecayMode->at(TauIndex),
								   my_genmatching_l2,
								   false  /// this is set to true for fake bakground
								   );
				   
				      //cout<<" sf : "<< applySf <<endl;
				      event_weight = event_weight * applySf;
				   
				      //makeTestPlot("f_dyll", 0,0,0,event_weight);
				      if( thirdLeptonVeto(MuIndex , TauIndex)  )
					{
					  nPassedThirdLepVeto_dyll+=event_weight;
					  //makeTestPlot("g_dyll", 0,0,0,event_weight);
					  //bool pbjv = (bJet_medium(EleIndex, TauIndex).size()==0) && (bJet_loose(EleIndex, TauIndex).size()<2);
					  if( pass_bjet_veto )
					    {
					      
					      nPassedBjetVeto_dyll+=event_weight;
					      //makeTestPlot("h_dyll", 0,0,0,event_weight);
					      double deltaR =  my_muP4.DeltaR(my_tauP4);
					      if(deltaR > 0.5 )
						{
						  nDeltaRPassed_dyll+=event_weight;
						  if(is_MC==false)event_weight=1.0;
						  //makeTestPlot("i_dyll", 0,0,0,event_weight);
						  if(debug)cout<<"this worked Line 374"<<endl;
						  fillHist("5_dyll",  MuIndex, TauIndex, false, event_weight);
						  fillUncPlots("5_dyll", MuIndex, TauIndex, false, event_weight);
						  double mT_muMet = TMass_F( my_muP4.Pt(), my_muP4.Phi(),
									     my_metP4.Pt(), my_metP4.Phi() );
						  if( mT_muMet < 50 )
						    {
						      fillHist("6_dyll", MuIndex, TauIndex, false, event_weight);
						      //makeTestPlot("j_dyll", 0,0,0,event_weight);
						    }
						  //higgs pt > 65,  visible mass < 125, met > 105
						  double mvis=(my_muP4+my_tauP4).M();
						  double higgsPt = pTvecsum_F(my_muP4, my_tauP4, my_metP4);
						  if(higgsPt > 65)
						    {
						      fillHist("7_dyll", MuIndex, TauIndex, false, event_weight);
						      if(mvis < 125)
							{
							  fillHist("8_dyll", MuIndex, TauIndex, false, event_weight);
							  if(my_metP4.Pt() > 105)
							    {
							      fillHist("9_dyll", MuIndex, TauIndex, false, event_weight);
							      // fillHist("9_dyll", MuIndex, TauIndex, false, event_weight);
							       fillUncPlots("9_dyll", MuIndex, TauIndex, false, event_weight);
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
	}


      ////// signal region -  isolated begin
      if(debug)cout<<"signal region -  isolated begin L523"<<endl;       
       
      Ztt_selector=false;
      ////// signal region -  isolated begin
      if(is_MC)
	event_weight=weight;
      else
	event_weight=1.0;
      //muCand.clear();  tauCand.clear();
       
      if(metFilters==0 )
	{
	  if(debug)cout<<"metfilters selected"<<endl;
	  if(is_MC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	  nMETFiltersPassed+=event_weight;
	  if(debug)cout<<"genweight applied"<<endl;
	  if(   passSingleTriggerPaths || passCrossTrigger )
	    {
	      nSingleTrgPassed+=event_weight;
	      if(debug)cout<<"trigger selected"<<endl;
	      //muCand = getMuCand(20,2.4);  ///// muons selected 
	      if( muCand.size() >0 ) 
		{ 
		  nGoodMuonPassed+=event_weight;
		  if(debug)cout<<"this worked Line 443"<<endl;
		   
		  //tauCand = getTauCand(30,2.3);
		  if( tauCand.size()>0 ) 
		    {
		      //nGoodTauPassed+=event_weight;
		      if(debug)cout<<"this worked Line 305"<<endl;
		      setMyEleTau(muCand[0], tauCand[0]);
		      //printTabSeparated("in Z->tt", MuIndex, TauIndex, my_muP4.Pt(), my_tauP4.Pt(), my_metP4.Pt(), my_njets );
		      // if( passDiMuonVeto(MuIndex)==true 
		      // 	  && eVetoZTTp001dxyz(MuIndex, TauIndex)
		      // 	  && mVetoZTTp001dxyz(MuIndex, TauIndex)
		      // 	  ) Ztt_selector=true;
		      // else Ztt_selector=false;
		       
		      if ( TriggerSelection(my_muP4, my_tauP4) )
			{
			  nGoodTauPassed+=event_weight;
			  //if(Ztt_selector) 
			  if(pass3rdLeptonVeto)
			    {
			      if ( muCharge->at(MuIndex) * tau_Charge->at(TauIndex) < 0  
				   &&  (if_DY_Genmatching(MuIndex, TauIndex)==1 || if_DY_Genmatching(MuIndex, TauIndex)==3)  )
				{
				  nGoodMuTauPassed+=event_weight;
				  //makeTestPlot("e", 0,0,0,event_weight);
				   
				  if ( MatchTriggerFilter(MuIndex, TauIndex) )
				    {
				       
				      if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
				       
				      //my_muP4 = my_muP4*muRC_sf;
				      applySf=1.0;
				      if(is_MC)
					applySf=  getScaleFactors( my_muP4.Pt(),
								   my_tauP4.Pt(),
								   my_muP4.Eta(),
								   my_tauP4.Eta(),
								   tau_DecayMode->at(TauIndex),
								   my_genmatching_l2,
								   false  /// this is set to true for fake bakground
								   );
				   
				      // if(debug)cout<<" sf : "<<getScaleFactors( EleIndex[0] , TauIndex[0] , false , is_MC , debug ) <<endl;
				      //if(my_muP4.Pt()<30)cout<<" sf : "<< applySf<< "  , mupt :"<<my_muP4.Pt() <<endl;
				      event_weight = event_weight * applySf;
				   
				      //makeTestPlot("f", 0,0,0,event_weight);
				      if( thirdLeptonVeto(MuIndex , TauIndex)  )
					{
					  nPassedThirdLepVeto+=event_weight;
					  //makeTestPlot("g", 0,0,0,event_weight);
					  //bool pbjv = (bJet_medium(EleIndex, TauIndex).size()==0) && (bJet_loose(EleIndex, TauIndex).size()<2);
					  if( pass_bjet_veto )
					    {
					      
					      nPassedBjetVeto+=event_weight;
					      //makeTestPlot("h", 0,0,0,event_weight);
					      double deltaR =  my_muP4.DeltaR(my_tauP4);
					      if(deltaR > 0.5 )
						{
						  nDeltaRPassed+=event_weight;
						  if(is_MC==false)event_weight=1.0;
						  //makeTestPlot("i", 0,0,0,event_weight);
						  if(debug)cout<<"this worked Line 374"<<endl;
						  fillHist("5",  MuIndex, TauIndex, false, event_weight);
						  fillUncPlots("5", MuIndex, TauIndex, false, event_weight);
						  double mT_muMet = TMass_F( my_muP4.Pt(), my_muP4.Phi(),
									     my_metP4.Pt(), my_metP4.Phi() );
						  if( mT_muMet < 50 )
						    {
						      fillHist("6", MuIndex, TauIndex, false, event_weight);
						      //makeTestPlot("j", 0,0,0,event_weight);
						    }
						  //higgs pt > 65,  visible mass < 125, met > 105
						  double mvis=(my_muP4+my_tauP4).M();
						  double higgsPt = pTvecsum_F(my_muP4, my_tauP4, my_metP4);
						  if(higgsPt > 65)
						    {
						      fillHist("7", MuIndex, TauIndex, false, event_weight);
						      if(mvis < 125)
							{
							  fillHist("8", MuIndex, TauIndex, false, event_weight);
							  if(my_metP4.Pt() > 105)
							    {
							      fillHist("9", MuIndex, TauIndex, false, event_weight);
							      // fillHist("9", MuIndex, TauIndex, false, event_weight);
							       fillUncPlots("9", MuIndex, TauIndex, false, event_weight);
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
	}
      /// signal region end
       
      ////// fake background region - antiisolated begin
      if(is_MC)
	event_weight=weight;
      else
	event_weight=1.0;
      //muCand.clear(); tauCand.clear();
      if(metFilters==0 )
	{
	  if(debug)cout<<"metfilters selected"<<endl;
	  if(is_MC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	  nMETFiltersPassed_fr+=event_weight;
	  if(debug)cout<<"genweight applied"<<endl;
	  if(   passSingleTriggerPaths || passCrossTrigger )
	    {
	      nSingleTrgPassed_fr+=event_weight;
	      if(debug)cout<<"trigger selected"<<endl;
	      //muCand = getMuCand(20,2.4);  ///// muons selected 
	      if( muCand.size() >0 ) 
		{ 
		  nGoodMuonPassed_fr+=event_weight;
		  if(debug)cout<<"this worked Line 443"<<endl;
		  
		  //tauCand = getAISRTauCand(30,2.3);
		  if( aisrtauCand.size()>0 ) 
		    {
		      //nGoodTauPassed_fr+=event_weight;
		      if(debug)cout<<"this worked Line 305"<<endl;
		      setMyEleTau(muCand[0], aisrtauCand[0]);
		      //printTabSeparated("in jet->tt", MuIndex, TauIndex, my_muP4.Pt(), my_tauP4.Pt(), my_metP4.Pt(), my_njets );
		      if ( TriggerSelection(my_muP4, my_tauP4) )
			{
			  double mt=TMass_F(my_muP4.Pt(),my_muP4.Phi()
					    ,my_metP4.Pt(), my_metP4.Phi());
			  double mvis=(my_muP4+my_tauP4).M();
			  double higgsPt = pTvecsum_F(my_muP4, my_tauP4, my_metP4);
			  double frac_tt=0.01; double frac_qcd=0.24; double frac_w=0.75; 
			  int category=eventCategory(MuIndex , TauIndex, higgsPt) ;
			  getFractions(category, mvis, frac_qcd, frac_w, frac_tt); /// this assigns right values for qcd, w and tt fractions
			  bool xtrg = false;
			  if( passCrossTrigger && my_muP4.Pt()<=25.0) xtrg=true;
			  else if ( my_muP4.Pt()>28.0) xtrg=false;
			  double newFF = FF_weights_withlpt.get_ff( my_tauP4.Pt(), mt, mvis
								    , 0 , my_muP4.Pt(), my_metP4.Pt()
								    , my_njets, xtrg
								    , frac_tt, frac_qcd, frac_w
								    , TString(" "));
				      
			  event_weight = event_weight*newFF; 
			  nGoodTauPassed_fr+=event_weight;
			  //if(Ztt_selector) 
			  if (pass3rdLeptonVeto)
			    {
			      if ( muCharge->at(MuIndex) * tau_Charge->at(TauIndex) < 0 )
				{
				  nGoodMuTauPassed_fr+=event_weight;
				  //makeTestPlot("e", 0,0,0,event_weight);
				  
				  if ( MatchTriggerFilter(MuIndex, TauIndex) )
				    {
				      
				      if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
				      
				      //my_muP4 = my_muP4*muRC_sf;
				      applySf=1.0;
				      if(is_MC)
					applySf=  getScaleFactors( my_muP4.Pt(),
								   my_tauP4.Pt(),
								   my_muP4.Eta(),
								   my_tauP4.Eta(),
								   tau_DecayMode->at(TauIndex),
								   my_genmatching_l2,
								   true  /// this is set to true for fake bakground
								   );
				      
				      // if(debug)cout<<" sf : "<<getScaleFactors( EleIndex[0] , TauIndex[0] , false , is_MC , debug ) <<endl;
				      //cout<<" sf : "<< applySf <<endl;
				      event_weight = event_weight * applySf;
				      //event_weight = event_weight* getFR(TauIndex);

				      //cout<<" newFF : "<< newFF <<endl;
				      //makeTestPlot("f_fr", 0,0,0,event_weight);
				      if( thirdLeptonVeto(MuIndex , TauIndex)  )
					{
					  nPassedThirdLepVeto_fr+=event_weight;
					  //makeTestPlot("g_fr", 0,0,0,event_weight);
					  //bool pbjv = (bJet_medium(EleIndex, TauIndex).size()==0) && (bJet_loose(EleIndex, TauIndex).size()<2);
					  if( pass_bjet_veto )
					    {
					      
					      nPassedBjetVeto_fr+=event_weight;
					      //makeTestPlot("h_fr", 0,0,0,event_weight);
					      double deltaR =  my_muP4.DeltaR(my_tauP4);
					      if(deltaR > 0.5 )
						{
						  nDeltaRPassed_fr+=event_weight;
						  //makeTestPlot("i_fr", 0,0,0,event_weight);
						  if(debug)cout<<"this worked Line 374"<<endl;
						  fillHist("5_fr",  MuIndex, TauIndex, true, event_weight);
						  fillUncPlots("5_fr", MuIndex, TauIndex, true, event_weight);						  
						  double mT_muMet = TMass_F( my_muP4.Pt(), my_muP4.Phi(),
									     my_metP4.Pt(), my_metP4.Phi() );
						  if( mT_muMet < 50 )
						    {
						      fillHist("6_fr", MuIndex, TauIndex, true, event_weight);
						      //makeTestPlot("j_fr", 0,0,0,event_weight);
						    }
						  //higgs pt > 65,  visible mass < 125, met > 105
						  double mvis=(my_muP4+my_tauP4).M();
						  double higgsPt = pTvecsum_F(my_muP4, my_tauP4, my_metP4);
						  if(higgsPt > 65)
						    {
						      fillHist("7_fr", MuIndex, TauIndex, true, event_weight);
						      if(mvis < 125)
							{
							  fillHist("8_fr", MuIndex, TauIndex, true, event_weight);
							  if(my_metP4.Pt() > 105)
							    {
							      fillHist("9_fr", MuIndex, TauIndex, true, event_weight);
							      // fillHist("9_fr", MuIndex, TauIndex, true, event_weight);
							       fillUncPlots("9_fr", MuIndex, TauIndex, true, event_weight);
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
   
   // fileName->cd();
   // map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
   // map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
   // for (; iMap1 != jMap1; ++iMap1)
   //   nplot1(iMap1->first)->Write();
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
       
       if( kinematic==true  &&  muonID==true &&  muonIso==true ){
	 tmpCand.push_back(iMu);
       }                                                                                      
     }                                                                                       
  return tmpCand;
  
}


std::vector<int> mutau_analyzer::getTauCand(double tauPtCut, double tauEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  TLorentzVector dau2;
  //Loop over taus      
  for(int iTau=0;iTau<nTau;iTau++)
    {
      dau2.SetPtEtaPhiE(tau_Pt->at(iTau),tau_Eta->at(iTau)
			,tau_Phi->at(iTau), tau_Energy->at(iTau)
			);
      if(is_MC)
	applyTauESCorrections(dau2, iTau, dau2);
      
      bool kinematic = false;
      bool trigger = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool mutau_separation=false;
      bool newDecayModeFinding=false;
      bool tau_reject=false;
      bool filter = false;
      if(  dau2.Pt() > tauPtCut 
	   && fabs( dau2.Eta() )< tauEtaCut
	   && tau_LeadChargedHadron_dz->at(iTau) < 0.2
	   )kinematic = true;
      if( tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1 ) tauIsolation=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byTightDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;

      if( kinematic==true    
	  && decayModeCut==true   
	  && tauIsolation==true 
	  && tau_reject==true   
	  && newDecayModeFinding==true
	  )
	{
	  tmpCand.push_back(iTau);
    	}                                                           
    }                                                                                       
  return tmpCand;
  
}
std::vector<int> mutau_analyzer::getAISRTauCand(double tauPtCut, double tauEtaCut){
  std::vector<int> tmpCand;  tmpCand.clear();
  TLorentzVector dau2;
  for(int iTau=0;iTau<nTau;iTau++) //Loop over taus
    {
      dau2.SetPtEtaPhiE(tau_Pt->at(iTau),tau_Eta->at(iTau)
			,tau_Phi->at(iTau), tau_Energy->at(iTau)
			);
      if(is_MC)
	applyTauESCorrections(dau2, iTau, dau2);

      bool kinematic = false;
      bool trigger = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool mutau_separation=false;
      bool newDecayModeFinding=false;
      bool tau_reject=false;
      bool filter = false;
      if( dau2.Pt() > tauPtCut 
	  && fabs( dau2.Eta())< tauEtaCut 
	  && tau_LeadChargedHadron_dz->at(iTau) < 0.2
  	  )kinematic = true;
      if( tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 && !(tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1) ) tauIsolation=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byTightDeepTau2017v2p1VSmu->at(iTau)==1 )tau_reject=true;
      
      
      if( kinematic==true    
	  && decayModeCut==true   
	  && tauIsolation==true 
	  && tau_reject==true   
	  && newDecayModeFinding==true
     	  )
	{
	  tmpCand.push_back(iTau);
    	}                                                           
    }                                                                                       
  return tmpCand;  
}
int mutau_analyzer::getZCand()
{
  if(!is_MC)
    return -1;
  for(int iMC=0;iMC<nMC;iMC++) //Loop over mc
    {
      if(fabs(mcPID->at(iMC))==23 && mcStatus->at(iMC)==62)
	return iMC;
    }
  return -1;
}
int mutau_analyzer::get_t_Cand()
{
  if(!is_MC)
    return -1;
  for(int iMC=0;iMC<nMC;iMC++) //Loop over mc
    {
      if(mcPID->at(iMC)==6 && mcStatus->at(iMC)==62)
	return iMC;
    }
  return -1;
}
int mutau_analyzer::get_tbar_Cand()
{
  if(!is_MC)
    return -1;
  for(int iMC=0;iMC<nMC;iMC++) //Loop over mc
    {
      if(mcPID->at(iMC)==-6 && mcStatus->at(iMC)==62)
        return iMC;
    }
  return -1;
}
std::vector<int> mutau_analyzer::getJetCand(int muIndex, int tauIndex){
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
         //&& (jetID->at(iJet)>>0&1)==1
         ) foundNoisyJets=true;

      if( jetPt->at(iJet) < 50 )
        {
          if(jetPUFullID->at(iJet)>>1&1==1 )
            passLoosePUID=true;
          else
            passLoosePUID=false;
        }
      else if (jetPt->at(iJet) > 50 )
        passLoosePUID=true;
      
      double lepton1Phi=muPhi->at(muIndex);
      double lepton1Eta= muEta->at(muIndex);
      double lepton2Phi=0;double lepton2Eta=0;
      lepton2Phi= tau_Phi->at(tauIndex); lepton2Eta=tau_Eta->at(tauIndex);
      double dr_jetEle=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton1Phi, lepton1Eta );
      double dr_jetTau=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton2Phi, lepton2Eta);
      if( dr_jetEle>0.5 && dr_jetTau>0.5 )
        drPassed=true;

      if(kinematic30 && !foundNoisyJets && passLoosePUID && drPassed==true)
        tmpCand.push_back(iJet);
    }
  return tmpCand;
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
bool mutau_analyzer::thirdLeptonVeto(int muIndex, int tauIndex){
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
  double deltaRm1=0; double deltaRm2=0; bool found_3rdmu=false;
  if(tmpCand.size() > 0 )
    { 
      deltaRm1 = delta_R(muPhi->at(muIndex),muEta->at(muIndex), elePhi->at(tmpCand[0]),  eleEta->at(tmpCand[0]));
      deltaRm2 = delta_R(tau_Phi->at(tauIndex),tau_Eta->at(tauIndex), elePhi->at(tmpCand[0]),  eleEta->at(tmpCand[0]));
      if(deltaRm1>0.5 && deltaRm2>0.5 ){
	return false;
      }
    }
  else
    return true;
  
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
vector<int> mutau_analyzer::bJet_medium(int muIndex, int tauIndex)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  double lepton1Phi=muPhi->at(muIndex);
  double lepton1Eta= muEta->at(muIndex);
  double lepton2Phi= tau_Phi->at(tauIndex); double lepton2Eta=tau_Eta->at(tauIndex);   
  double dr_jetEle=0.0; double dr_jetTau=0.0; 
  
  for(int iJets=0; iJets<nJet ; iJets++){
    bool particles_separated=false;
    bool kinematic = false;
    bool passLoosePUID=false;
    dr_jetEle =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton1Phi, lepton1Eta );
    dr_jetTau =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton2Phi, lepton2Eta);
    if( dr_jetEle>0.5 && dr_jetTau>0.5)
      particles_separated=true;
    if( jetPt->at(iJets) > 25
	&& abs(jetEta->at(iJets)) < 2.4 
	&& jetID->at(iJets)>>0&1==1 
       	&& (jetDeepCSVTags_b->at(iJets) + jetDeepCSVTags_bb->at(iJets)) > 0.4941
	)
      kinematic=true;
    if( jetPt->at(iJets)<50 )
      {
	if(jetPUFullID->at(iJets)>>1&1==1 )
	  passLoosePUID=true;
	else
	  passLoosePUID=false;
      }
    else if (jetPt->at(iJets) > 50 )
      passLoosePUID=true;
    
    if(particles_separated && kinematic )
      tmpJetCand.push_back(iJets);
  }
  return tmpJetCand;
}

vector<int> mutau_analyzer::bJet_loose(int muIndex, int tauIndex)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  double lepton1Phi=muPhi->at(muIndex);
  double lepton1Eta= muEta->at(muIndex);
  double lepton2Phi= tau_Phi->at(tauIndex); double lepton2Eta=tau_Eta->at(tauIndex);   
  double dr_jetEle=0.0; double dr_jetTau=0.0; 
  
  for(int iJets=0; iJets<nJet ; iJets++){
    bool kinematic = false;
    bool passLoosePUID=false;
    bool particles_separated=false;
    dr_jetEle =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton1Phi, lepton1Eta );
    dr_jetTau =delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton2Phi, lepton2Eta);
    if( dr_jetEle>0.5 && dr_jetTau>0.5)
      particles_separated=true;
      
    if( jetPt->at(iJets) > 25
	&& abs(jetEta->at(iJets)) < 2.4 
	&& jetID->at(iJets)>>0&1==1 
	&& (jetDeepCSVTags_b->at(iJets) + jetDeepCSVTags_bb->at(iJets)) > 0.1522 
	)
      kinematic=true;
    if( jetPt->at(iJets)<50 )
      {
        if(jetPUFullID->at(iJets)>>1&1==1 )
          passLoosePUID=true;
        else
          passLoosePUID=false;
      }
    else if (jetPt->at(iJets) > 50 )
      passLoosePUID=true;
    
    if(particles_separated && kinematic )
      tmpJetCand.push_back(iJets);
    
  }
  return tmpJetCand;
}

bool passBjetVeto_medium(int eleIndex, int tauIndex){
  return true;
}
bool passBjetVeto_loose(int eleIndex, int tauIndex){
  return true;
}

bool passBjetVeto(int eleIndex, int tauIndex){
  return passBjetVeto_medium(eleIndex,tauIndex) && passBjetVeto_loose(eleIndex,tauIndex);
  
}
bool mutau_analyzer::passDiMuonVeto(int muIndex)
{
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool veto = true;
  bool awayFromEverything=true;
  for(int iMu=0;iMu<nMu;iMu++)
    {
      bool kinematic = false;
      if( muPt->at(iMu) > 15
	  && fabs(muEta->at(iMu))< 2.4
	  && fabs(muD0->at(iMu)) < 0.045
	  && fabs(muDz->at(iMu)) < 0.2
       	  ) kinematic = true;
      bool muonId =false;
      if( muIDbit->at(iMu)>>6&1==1) muonId =true; // pf muon iso veryloose
      bool relative_iso = false;    
      float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
      if( relMuIso < 0.3 ) relative_iso = true;
      if( kinematic && muonId && relative_iso ){
	tmpCand.push_back(iMu);
      }	
    }
  std::vector<int> iMuPlus;  iMuPlus.clear(); 
  std::vector<int> iMuMinus; iMuMinus.clear();
  for(int i=0; i<tmpCand.size(); i++){
    if(muCharge->at(tmpCand[i]) < 0) iMuMinus.push_back(tmpCand[i]);
    if(muCharge->at(tmpCand[i]) > 0) iMuPlus.push_back(tmpCand[i]);
  }
  if(iMuPlus.size()>0 && iMuMinus.size()>0){
    double deltaR= delta_R(muPhi->at(iMuMinus[0]), muEta->at(iMuMinus[0])
			   , muPhi->at(iMuPlus[0]), muEta->at(iMuPlus[0]));
    if (deltaR > 0.15 && muCharge->at(iMuPlus[0])*muCharge->at(iMuMinus[0])<0) {
      return false;
    }
  }
  return true;
  
}

bool mutau_analyzer::eVetoZTTp001dxyz(int muIndex, int tauIndex){
  std::vector<int> tmpCand;  tmpCand.clear();
  std::vector<int> output;  output.clear();
  bool awayFromEverything = true;   int tmpEleIndex=-1;
  //Loop over electrons      
  for(int iEle=0;iEle<nEle;iEle++)
    {
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
	  double deltaR_ee = delta_R(muPhi->at(muIndex), muEta->at(muIndex), elePhi->at(tmpCand[i]), eleEta->at(tmpCand[i]));
	  if(! (deltaR_et>0.0001 && deltaR_ee>0.0001))
	    output.push_back(i);
	}
    }
  if(output.size() >0 )
    return false;
  else
    return true;
    
}
bool mutau_analyzer::mVetoZTTp001dxyz(int muIndex, int tauIndex){
  std::vector<int> tmpCand; tmpCand.clear();
  std::vector<int> output;  output.clear();
  bool awayFromEverything = true;   int tmpMuIndex=-1;
  //Loop over muons
  for(int iMu=0; iMu < nMu;iMu++)
    {
      //if (iMu==muIndex)continue;
      bool kinematic = false;
      if( (*muPt)[iMu] > 10.0  && fabs((*muEta)[iMu])< 2.4 && (*muD0)[iMu] < 0.045 && (*muDz)[iMu] < 0.2 ) kinematic = true;
      bool muonId =false;
      if( muIDbit->at(iMu)>>8&1==1) muonId =true;
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
      for(int i=0;i<tmpCand.size();i++)
        {
	  deltaRm1 = delta_R(muPhi->at(muIndex),muEta->at(muIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
	  deltaRm2 = delta_R(tau_Phi->at(tauIndex),tau_Eta->at(tauIndex), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
	  if(! (deltaRm1>0.0001 && deltaRm2>0.0001) )
	    output.push_back(i);
	}
    }
  if(output.size() >1 )
    return false;
  else
    return true;
  
}
int mutau_analyzer::myGenMaching(int tauIndex)
{
  if(is_MC==false)
    return 0;

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
int mutau_analyzer::myGenMaching1(int muIndex)
{
  if(is_MC==false)
    return 0;
  double recotau_eta=muEta->at(muIndex);
  double recotau_phi=muPhi->at(muIndex);
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



int mutau_analyzer::getGenMu(){
  std::vector<int> tmpCand;
  tmpCand.clear();
  //int count1=0; int count2=0;
  for(int imc=0; imc<nMC; imc++){
    if( (mcStatus->at(imc)==23 || mcStatus->at(imc)==33) && fabs(mcPID->at(imc))==13 )
      tmpCand.push_back(imc);
  }
  if (tmpCand.size()>0)
    {
      //cout<<"gen mu found "<<tmpCand[0]<<endl;
      return tmpCand[0];
    }
  return -99;
}

float mutau_analyzer::exponential(float x,float a,float b,float c) {
  return a * TMath::Exp(-b * x) + c;
}
double mutau_analyzer::exponential( double pT) {
  //return TMath::Exp(0.088 - 0.00087*pT + 0.00000092*pow(pT,2) ) ; 
  return  TMath::Exp(0.0615-0.0005*pT );
}

double mutau_analyzer::getScaleFactors(  double mupt, double taupt, double mueta, double taueta, int taudm, int tauGenMatch, bool isFakebkg)
{
  double rv_sf=1.0;
  double sf_IsoEff = 1.0; 
  double sf_muTrg = 1.0;
  double sf_muID = 1.0;
  double sf_tauidSF_m = 1.0;
  double sf_tauTrg = 1.0;
  double sf_tauidSF_vvvl = 1.0;
  double sf_tauesSF = 1.0;
  double sf_fakeEle = 1.0; double sf_fakeMu = 1.0;
  double sf_taufesSF = 1.0;
  
  double recoMuonPt=0.0;
  if ( mupt < 120)
    recoMuonPt=mupt;
  else
    recoMuonPt = 119;
  sf_muID = h_muIDSF->GetBinContent(h_muIDSF->GetXaxis()->FindBin(recoMuonPt),h_muIDSF->GetYaxis()->FindBin(abs(mueta)));
  sf_IsoEff = h_muIsoSF->GetBinContent(h_muIsoSF->GetXaxis()->FindBin(recoMuonPt),h_muIsoSF->GetYaxis()->FindBin(abs(mueta)));  
  sf_muTrg = h_muTrgSF->GetBinContent(h_muTrgSF->GetXaxis()->FindBin(mupt),h_muTrgSF->GetYaxis()->FindBin(abs(mueta)));

  if(  tauGenMatch>=5 )
    {
      sf_tauidSF_m = h_tauidSF_m->GetBinContent(h_tauidSF_m->GetXaxis()->FindBin(taudm));
      //sf_tauidSF_m = fn_tauIDSF_m->Eval(tau_Pt->at(TauIndex));
      sf_tauidSF_vvvl = h_tauidSF_vvvl->GetBinContent(h_tauidSF_vvvl->GetXaxis()->FindBin(taudm));
      //sf_tauidSF_vvvl = fn_tauIDSF_vvl->Eval(tau_Pt->at(TauIndex));
      sf_tauesSF = h_tauesSF->GetBinContent(h_tauesSF->GetXaxis()->FindBin(taudm));
    }

  // if(tauGenMatch==1 ||tauGenMatch==3 ) /// electrons to pass Deep tau
  //   {
  //     if(taudm==0)
  // 	{
  // 	  if(abs(taueta) < 1.479 ) sf_fakeEle=0.98;
  // 	  if(abs(taueta) > 1.479 ) sf_fakeEle=0.80;
  // 	}
  //     if(taudm==1)
  // 	{
  // 	  if(abs(taueta) < 1.479 ) sf_fakeEle=1.07;
  // 	  if(abs(taueta) > 1.479 ) sf_fakeEle=0.64;
  // 	}
  //   }
  // if(tauGenMatch==2 ||tauGenMatch==4){  ///  muons to pass deep tau 
  //   if(taudm==0)
  //     {
  // 	if(abs(taueta) < 0.4 ) sf_fakeMu=1.17;
  // 	if(abs(taueta) > 0.4 && abs(taueta) < 0.8 ) sf_fakeMu=1.29;
  // 	if(abs(taueta) > 0.8 && abs(taueta) < 1.2 ) sf_fakeMu=1.14;
  // 	if(abs(taueta) > 1.2 && abs(taueta) < 1.7 ) sf_fakeMu=0.93;
  // 	if(abs(taueta) > 1.7 && abs(taueta) < 2.3 ) sf_fakeMu=1.61;
  //     }
  //   if(taudm==1)
  //     {
  // 	if(abs(taueta) > 0.0 && abs(taueta) < 0.4 ) sf_fakeMu=0.69;
  //     }
  // }
  if(tauGenMatch==2 ||tauGenMatch==4){  ///  muons to pass deep tau
    if(abs(taueta) < 0.4 ) sf_fakeMu=1.17;
    else if(abs(taueta) > 0.4 && abs(taueta) < 0.8 ) sf_fakeMu=1.29;
    else if(abs(taueta) > 0.8 && abs(taueta) < 1.2 ) sf_fakeMu=1.14;
    else if(abs(taueta) > 1.2 && abs(taueta) < 1.7 ) sf_fakeMu=0.93;
    else if(abs(taueta) > 1.7 && abs(taueta) < 2.3 ) sf_fakeMu=1.61;
  }
  if(tauGenMatch==1 ||tauGenMatch==3){  ///  muons to pass deep tau
    if(abs(taueta) < 1.460 ) sf_fakeMu=1.09;
    else if  (abs(taueta) >=1.559 ) sf_fakeMu=1.19;
  }
  double higgsPt = pTvecsum_F(my_muP4, my_tauP4, my_muP4);
  double higgPt_weight=1.0;
  if (my_njets==0)
    higgPt_weight = gr_NNLOPSratio_pt_mcatnlo_0jet->Eval(min(higgsPt,125.0));
  else if (my_njets==1)
    higgPt_weight = gr_NNLOPSratio_pt_mcatnlo_1jet->Eval(min(higgsPt,625.0));
  else if (my_njets==2)
    higgPt_weight = gr_NNLOPSratio_pt_mcatnlo_2jet->Eval(min(higgsPt,800.0));
  else if (my_njets>=3)
    higgPt_weight = gr_NNLOPSratio_pt_mcatnlo_3jet->Eval(min(higgsPt,925.0));
  else
    higgPt_weight = 1.0;
  //  if (higgPt_weight>1.0) higgPt_weight = 1.0;
  
  double tauPtCheck=taupt;
  if(taupt > 450 ) tauPtCheck = 450;
  else if ( taupt < 20 )  tauPtCheck = 20;
  
  if(taudm==0)  sf_tauTrg= h_tauTrgSF_dm0->GetBinContent(h_tauTrgSF_dm0->GetXaxis()->FindBin(tauPtCheck));
  if(taudm==1)  sf_tauTrg= h_tauTrgSF_dm1->GetBinContent(h_tauTrgSF_dm1->GetXaxis()->FindBin(tauPtCheck));
  if(taudm==10) sf_tauTrg= h_tauTrgSF_dm10->GetBinContent(h_tauTrgSF_dm10->GetXaxis()->FindBin(tauPtCheck));
  if(taudm==11) sf_tauTrg= h_tauTrgSF_dm11->GetBinContent(h_tauTrgSF_dm11->GetXaxis()->FindBin(tauPtCheck));
  
  ///// btag efficiency
  double weight_btagSF = btag_sf;
  //cout<<"btag_sf : "<<weight_btagSF<<endl;

   //// ele trigger, id scale factors from RooWorkspace
  w->var("m_pt")->setVal(mupt);
  w->var("m_eta")->setVal(mueta);
  w->var("t_pt")->setVal(taupt);
  w->var("t_mvadm")->setVal(taudm);
  double m_tracking = w->function("m_trk_ratio")->getVal();
  double m_IDiso = w->function("m_idiso_ic_ratio")->getVal();
  double mu_singletrg = w->function("m_trg_ic_ratio")->getVal();
  double mu_crosstrg = w->function("m_trg_20_ic_ratio")->getVal();
  double tau_crossttrg = w->function("t_trg_ic_deeptau_medium_mvadm_mutau_ratio")->getVal();
  double t_deepid = w->function("t_deeptauid_mvadm_medium")->getVal();
  double zptmass_weight = 1.0;
  if(found_DYjet_sample)
    zptmass_weight= get_zptmass_weight();
  if(zptmass_weight>1) zptmass_weight = 1.0;

  double top_pt_weight=1.0;
  if(found_TTbar_sample){
    int t_index = get_t_Cand(); int tbar_index = get_tbar_Cand();
    if( t_index >-1 && tbar_index > -1 ){
      double pttop1= mcPt->at(t_index); double pttop2= mcPt->at(tbar_index);
      if (pttop1>400) pttop1 = 400;
      if (pttop2>400) pttop2 = 400;
      top_pt_weight = sqrt( exponential(pttop1) * exponential(pttop2) );
      //cout<<"top_pt_weight = "<<top_pt_weight<<endl;
    } 
  }
  if ( mupt<25.0 ) { mu_singletrg=1.0; }
  else  { mu_crosstrg=1.0; tau_crossttrg = 1.0; sf_tauTrg = 1.0; tau_crossttrg=1.0; }

  double sf_htt_workspace=  m_tracking * m_IDiso *  mu_singletrg * mu_crosstrg * tau_crossttrg * t_deepid;// * zptmass_weight;

  rv_sf = sf_tauidSF_m * sf_tauTrg * sf_fakeMu * sf_htt_workspace * top_pt_weight * higgPt_weight ;
  rv_sf = rv_sf * MuonIDIso.get_ScaleFactor(my_muP4.Pt(), my_muP4.Eta());
  if (  my_muP4.Pt()<25.0 )
    {
      rv_sf = rv_sf* CrossTriggerSF.get_ScaleFactor(my_muP4.Pt(), my_muP4.Eta()) ;
      //rv_sf = rv_sf * MuonIDIso.get_ScaleFactor(my_muP4.Pt(), my_muP4.Eta());
    }
  else 
    rv_sf = rv_sf*IsoMu24or27SF.get_ScaleFactor(my_muP4.Pt(), my_muP4.Eta());
  // if(isFakebkg)
  //   rv_sf = rv_sf * sf_tauidSF_vvvl;
  //cout<<" "<< MuonIDIso.get_ScaleFactor(my_muP4.Pt(), my_muP4.Eta()) <<" \t" <<  IsoMu24or27SF.get_ScaleFactor(my_muP4.Pt(), my_muP4.Eta()) <<endl;  

  return rv_sf;
}

bool mutau_analyzer::TriggerSelection(TLorentzVector muP4, TLorentzVector tauP4){

  if( muP4.Pt() > 21.0 && muP4.Pt() < 25.0 && tauP4.Pt()>32.0  && fabs(muP4.Eta())< 2.1  && fabs(tauP4.Eta())< 2.1 ){
    if( HLTTau>>0&1 == 1)
      return true;
    else
      return false;
  }
  else if ( muP4.Pt() > 25.0  && fabs(muP4.Eta())< 2.1  && tauP4.Pt()>30.0 && fabs(tauP4.Eta())< 2.3 ){
    if(  ( HLTEleMuX>>21&1 == 1 && muP4.Pt() > 28.0  ) //HLT_IsoMu27_v
	 || ( HLTEleMuX>>60&1 == 1 && muP4.Pt() > 25.0  ) //HLT_IsoMu24_v
	 )
      return true;
  }
  else
    return false;
  
}

bool mutau_analyzer::MatchTriggerFilter(int muIndex, int tauIndex)
{
  // std::vector<int> tmpJetCand;
  // tmpJetCand.clear();
  // bool passFilter = true;
  // bool muTriggerFilterMatch=false;
  // int nMuTriggerFilterMatch=0;
  // bool tauTriggerFilterMatch=false;
  // int nTauTriggerFilterMatch=0;
  // for(int ifilter=33;ifilter<56;ifilter++)
  //   {
  //     if(muFiredTrgs->at(muIndex)>>ifilter&1==1)
  // 	{
  // 	  muTriggerFilterMatch=true;
  // 	  nMuTriggerFilterMatch++;
  // 	}
  //   }
  // for(int ifilter=0;ifilter<18;ifilter++)
  //   {
  //     if(tauFiredTrgs->at(tauIndex)>>ifilter&1==1)
  // 	{
  // 	  tauTriggerFilterMatch=true;
  // 	  nTauTriggerFilterMatch++;
  // 	}
  //   }
  // if( (HLTEleMuX>>21&1 == 1 && nMuTriggerFilterMatch==23 )
  //     || (HLTEleMuX>>60&1 == 1 && nMuTriggerFilterMatch==23 ) 
  //     || (HLTTau>>0&1 ==1  && nMuTriggerFilterMatch==23 )
  //     )
  //   passFilter=true;
    
  // return passFilter;
  return true;

}



double  mutau_analyzer::getFR(int tauIndex){
  double frWeight=1.0;
  // double tau_FR = 1.0;
  // double tauPt=0.0;
  // if( tau_Pt->at(tauIndex) < 120 )
  //   tauPt=tau_Pt->at(tauIndex);
  // else
  //   tauPt=119.0;
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
  return frWeight;
}
void mutau_analyzer::fillUncPlots( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  //////////////////////  fill unccertainity plots ////////////////////////
  string hNumber = histNumber;

  ///////// jet-tau fakes unc
  if(isFakeBkg){
    jetFakeUnc( histNumber, muIndex, tauIndex, isFakeBkg, event_weight);
  }
  if(!is_MC) return;
  if(isFakeBkg) return;
  /////// tauh ID unc
  tauIDunc( histNumber, muIndex, tauIndex, true, event_weight);

  ////////// tau es unc
  tauESunc( histNumber, muIndex, tauIndex, true, event_weight);

  /////// trigger unc
  tauTRGunc( histNumber, muIndex, tauIndex, true, event_weight);
  leptonTRGunc( histNumber, muIndex, tauIndex, true, event_weight);
  /////// prefiring unc
  // prefiringweight , prefiringweightup, prefiringweightdown // only for 2017
  fillHist(histNumber+"_CMS_Prefiring_up", MuIndex, TauIndex, true, event_weight*prefiringweightup/prefiringweight);
  fillHist(histNumber+"_CMS_Prefiring_down", MuIndex, TauIndex, true, event_weight*prefiringweightdown/prefiringweight);

  /////// antilepton

  /////// muon -> tau es
  muonMissID(histNumber, MuIndex, TauIndex, true, event_weight);

  /////// muon es
  muonES(histNumber, MuIndex, TauIndex, true, event_weight);
  
  ////// JES
  JESUnc(histNumber, MuIndex, TauIndex, true, event_weight);

  ////// JER
  JERUnc(histNumber, MuIndex, TauIndex, true, event_weight);
  
  //// dy_shape
  dyShape(histNumber, MuIndex, TauIndex, true, event_weight);
  //// ttbar shape
  ttbarShape(histNumber, MuIndex, TauIndex, true, event_weight);
  
  ////// met sys
  metSysUnc(histNumber, MuIndex, TauIndex, true, event_weight);
  metClusteredUnc(histNumber, MuIndex, TauIndex, true, event_weight);

    
}

void mutau_analyzer::jetFakeUnc( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
 
  string jetFakesUnc[15] = {"rawFF_mt_qcd_unc1_2017", "rawFF_mt_qcd_unc2_2017", "rawFF_mt_tt_unc1_2017", "rawFF_mt_tt_unc2_2017", "rawFF_mt_w_unc1_2017", "rawFF_mt_w_unc2_2017", "FF_closure_OSSS_mvis_mt_qcd_2017", "FF_closure_lpt_mt_qcd", "FF_closure_lpt_mt_tt", "FF_closure_lpt_mt_w", "FF_closure_lpt_xtrg_mt_qcd_2017", "FF_closure_lpt_xtrg_mt_tt_2017", "FF_closure_lpt_xtrg_mt_w_2017", "FF_closure_mt_mt_w_unc1_2017", "FF_closure_mt_mt_w_unc2_2017"};

  string uncs[2] = {"up", "down"};

  double mt=TMass_F(my_muP4.Pt(),my_muP4.Phi()
		    ,my_metP4.Pt(), my_metP4.Phi());
  double mvis=(my_muP4+my_tauP4).M();
  double higgsPt = pTvecsum_F(my_muP4, my_tauP4, my_metP4);
  double frac_tt=0.01; double frac_qcd=0.24; double frac_w=0.75; 
  int category=eventCategory(MuIndex , TauIndex, higgsPt) ;
  getFractions(category, mvis, frac_qcd, frac_w, frac_tt); /// this assigns right values for qcd, w and tt fraction
  bool xtrg = false;
  if( (HLTTau>>0&1 == 1 || HLTTau>>14&1 == 1) && my_muP4.Pt()<=25.0) xtrg=true;
  else if ( my_muP4.Pt()>25.0) xtrg=false;

  double FF_nominal = FF_weights_withlpt.get_ff( my_tauP4.Pt(), mt, mvis
						 , 0 , my_muP4.Pt(), my_muP4.Pt()
						 , my_njets, xtrg
						 , frac_tt, frac_qcd, frac_w
						 ,TString(" "));

  //cout<<" \t "<< " FF_nominal = "<< FF_nominal <<endl;
  for(int u = 0; u < 2; u++) {
    for(int i = 0; i < 15; i++) {
      //cout << jetFakesUnc[i] << "\n";
      string uncName = jetFakesUnc[i] + "_" + uncs[u];
      
      double FF_unc = FF_weights_withlpt.get_ff( my_tauP4.Pt(), mt, mvis
						 , 0 , my_muP4.Pt(), my_muP4.Pt()
						 , my_njets, xtrg
						 , frac_tt, frac_qcd, frac_w
						 ,TString(uncName));
      double weight = event_weight * FF_unc / FF_nominal;
      if ( event_weight==0 || FF_nominal==0|| FF_unc==0)
	weight = 0;
      fillHist(histNumber+"_CMS_"+uncName, MuIndex, TauIndex, true, weight);
    }
  }
}
void mutau_analyzer::JESUnc( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  // up =  1 + jetJECUnc  , down = 1 - jetJECUnc
  vector<float> tmp_jetPt;
  float tmp_pfMET = pfMET;
  float tmp_pfMETPhi = pfMETPhi;
  for (float pt : (*jetPt)) tmp_jetPt.push_back(pt);

  for (int i = 0; i < nJet; i++)
    jetPt->at(i) = tmp_jetPt[i]*(1 + 1*jetJECUnc->at(i));

  pfMET = pfMET_T1JESUp;
  pfMETPhi = pfMETPhi_T1JESUp;
  setMyEleTau( muIndex, tauIndex);
  fillHist(histNumber+"_CMS_JES_up", MuIndex, TauIndex, true, event_weight);
  
  for (int i = 0; i < nJet; i++)
    jetPt->at(i) = tmp_jetPt[i]*(1 - 1*jetJECUnc->at(i));
  pfMET = pfMET_T1JESDo;
  pfMETPhi = pfMETPhi_T1JESDo;
  setMyEleTau( muIndex, tauIndex);
  fillHist(histNumber+"_CMS_JES_down", MuIndex, TauIndex, true, event_weight);
  
  /////////// we need to reset the jet, met variable
  for (int i = 0; i < nJet; i++)
    jetPt->at(i) = tmp_jetPt[i];
  pfMET = tmp_pfMET;
  pfMETPhi = tmp_pfMETPhi;
      
}
void mutau_analyzer::JERUnc( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  
  // if (isData || isSignal) return; ////// why
  vector<float> tmp_jetPt;
  float tmp_pfMET = pfMET;
  float tmp_pfMETPhi = pfMETPhi;
  for (float pt : (*jetPt)) tmp_jetPt.push_back(pt);

  for (int i = 0; i < nJet; i++)
    jetPt->at(i) = jetRawPt->at(i) * jetP4SmearUp->at(i);
    
  pfMET = pfMET_T1JERUp;
  setMyEleTau( muIndex, tauIndex);
  fillHist(histNumber+"_CMS_JER_up", MuIndex, TauIndex, true, event_weight);
  
  for (int i = 0; i < nJet; i++)
    jetPt->at(i) = jetRawPt->at(i) * jetP4SmearDo->at(i);
  pfMET = pfMET_T1JERDo;
  setMyEleTau( muIndex, tauIndex);
  fillHist(histNumber+"_CMS_JER_down", MuIndex, TauIndex, true, event_weight);
  

  /////////// we need to reset the jet, met variable
  for (int i = 0; i < nJet; i++)
    jetPt->at(i) = tmp_jetPt[i];
  pfMET = tmp_pfMET;
  pfMETPhi = tmp_pfMETPhi;

}
void mutau_analyzer::tauESunc( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  float pt_low  = 34;
  float pt_high = 170;
  int taudm = tau_DecayMode->at(tauIndex);
  string uncName = ""; string dmName = "";

  float tes = get_BinContent(h_tauesSF, taudm); float err = 0;
  if(my_tauP4.Pt() >= pt_high){
    err = get_BinError(h_tauesSF_highpt, taudm);
  }
  else if (my_tauP4.Pt() >pt_low ){
    float err_high = get_BinError(h_tauesSF_highpt, taudm);
    float err_low  = get_BinError(h_tauesSF, taudm);
    err = err_low + (err_high-err_low)/(pt_high - pt_low)*(my_tauP4.Pt() - pt_low);
  }
  else {
    err =get_BinError(h_tauesSF, taudm);
  }
  
  for(int i=0; i<4; i++){
    dmName = tauDM_map[i].first;
    string histName = histNumber+"_CMS_scale_t_"+dmName+"_2017";
    if(tauDM_map[i].second != taudm) err = 0;
    //cout<<"filling for "<<histName<<" dm = "<<tauDM_map[i].second<<endl;
    float weight = event_weight * (tes + err) / tes;
    fillHist(histName+"_up", MuIndex, TauIndex, true, weight);
    weight = event_weight * (tes - err) / tes;
    fillHist(histName+"_down", MuIndex, TauIndex, true, weight);
  }
}

void mutau_analyzer::tauTRGunc( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  int taudm = tau_DecayMode->at(tauIndex);
  double tauPtCheck= min(450.0, max(20.0, my_tauP4.Pt()));
  string uncName = "";
  float trg = 1.0; float err = 0;
  if(taudm==0){
    trg = get_BinContent(h_tauTrgSF_dm0, tauPtCheck);
    err = get_BinError(h_tauTrgSF_dm0, tauPtCheck);
  }
  else if(taudm==1){
    trg = get_BinContent(h_tauTrgSF_dm1, tauPtCheck);
    err = get_BinError(h_tauTrgSF_dm1, tauPtCheck);
  }
  else if(taudm==10){
    trg = get_BinContent(h_tauTrgSF_dm10, tauPtCheck);
    err = get_BinError(h_tauTrgSF_dm10, tauPtCheck);
  }
  else{
    trg = get_BinContent(h_tauTrgSF_dm11, tauPtCheck);
    err = get_BinError(h_tauTrgSF_dm11, tauPtCheck);
  }
  for(int i=0; i<4; i++){
    string histName = histNumber+"_CMS_scale_tautrg_"+tauDM_map[i].first+"_2017";
    if(tauDM_map[i].second != taudm) err = 0;
    float weight = event_weight * (trg + err) / trg;
    fillHist(histName+"_up", MuIndex, TauIndex, true, weight);
    weight = event_weight * (trg - err) / trg;
    fillHist(histName+"_down", MuIndex, TauIndex, true, weight);
  }
}

void mutau_analyzer::leptonTRGunc( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){

  string uncName = "";
  string mutrg[2] = {"CMS_mutautrg_2017", "CMS_singlemutrg_2017"};
  float weight = 1.0;
  w->var("m_pt")->setVal(my_muP4.Pt());
  w->var("m_eta")->setVal(my_muP4.Eta());
  w->var("t_pt")->setVal(my_tauP4.Pt());
  w->var("t_mvadm")->setVal(tau_DecayMode->at(tauIndex));
  double m_tracking = w->function("m_trk_ratio")->getVal();
  double m_IDiso = w->function("m_idiso_ic_ratio")->getVal();
  double mu_singletrg = w->function("m_trg_ic_ratio")->getVal();
  double mu_crosstrg = w->function("m_trg_20_ic_ratio")->getVal();
  
  /// fill cross trg
    {
      uncName = "CMS_mutautrg_2017";
      double err = 0.0; 
      double trg = 1.0; 
      if( my_muP4.Pt()<25.0 && my_tauP4.Pt()>32.0)
	{ trg = mu_crosstrg; err = 0.02;}
      weight = event_weight * (trg * (1+err)) / trg;
      fillHist(histNumber+"_"+uncName+"_up", MuIndex, TauIndex, true, weight);
      weight = event_weight * (trg * (1-err)) / trg;
      fillHist(histNumber+"_"+uncName+"_down", MuIndex, TauIndex, true, weight);
    }
    {
      uncName = "CMS_singlemutrg_2017";
      double err = 0.0;
      double trg = 1.0;
      if( my_muP4.Pt() > 25.0)
	{ trg = mu_singletrg; err = 0.02; }
      weight = event_weight * (trg * (1+err)) / trg;
      fillHist(histNumber+"_"+uncName+"_up", MuIndex, TauIndex, true, weight);
      weight = event_weight * (trg * (1-err)) / trg;
      fillHist(histNumber+"_"+uncName+"_down", MuIndex, TauIndex, true, weight);
    }
}

void mutau_analyzer::electronMissID( string histNumber , int eleIndex, int tauIndex, bool isFakeBkg, float event_weight){

  if( !(my_genmatching_l2==1 || my_genmatching_l2==3 || my_genmatching_l1==1 || my_genmatching_l1==3) )
    return;

  string regions[3] = {"barrel" , "endcap"};
  string taudm = "";
  float err = 0.01;
  double weight = 1.0;
  double sf_fes =  getTauFES(tauIndex, "");
  double sf_fes_up =  getTauFES(tauIndex, "up");
  double sf_fes_down =  getTauFES(tauIndex, "down");
  
  for(int i=0; i<2; i++){
    taudm = tauDM_map[i].first;
    if(tauDM_map[i].second != tau_DecayMode->at(tauIndex)) err = 0;
      if (abs(my_tauP4.Eta())< 1.5) 
    { //////// region = barrel
      string histName = histNumber + "_CMS_scale_efaket_" + taudm + "_" + "barrel" +"_2017";
      weight = event_weight * sf_fes_up / sf_fes;
      fillHist(histName+"_up", MuIndex, TauIndex, true, weight);
      weight = event_weight * sf_fes_down / sf_fes;
      fillHist(histName+"_down", MuIndex, TauIndex, true, weight);
    }
  else
    {
      string histName = histNumber + "_CMS_scale_efaket_" + taudm + "_" + "endcap" +"_2017";
      weight = event_weight * sf_fes_up / sf_fes;
      fillHist(histName+"_up", MuIndex, TauIndex, true, weight);
      weight = event_weight * sf_fes_down / sf_fes;
      fillHist(histName+"_down", MuIndex, TauIndex, true, weight);
    }
  }

}
void mutau_analyzer::electronES( string histNumber , int eleIndex, int tauIndex, bool isFakeBkg, float event_weight){

  string histName = histNumber + "_CMS_scale_e_Scale_2017";  
  float err = eleEnergySigma->at(eleIndex);
  float weight  = 1.0;

  weight = event_weight * (1 + err);
  fillHist(histName+"_up", MuIndex, TauIndex, true, weight);
  weight = event_weight * (1 - err);
  fillHist(histName+"_down", MuIndex, TauIndex, true, weight);

}
void mutau_analyzer::muonES( string histNumber , int eleIndex, int tauIndex, bool isFakeBkg, float event_weight){

  string regions[3] = {"etalt1p2", "eta1p2to2p1", "etam2p1to2p4" };
  bool applySys[3] = {0, 0, 0};
  float err = 0;
  double weight = 1.0;

  if (abs(my_tauP4.Eta()) < 1.2)
    {
      err = 0.004;
      applySys[0] = 1;
    }
  else if (abs(my_tauP4.Eta()) > 1.2  && abs(my_tauP4.Eta()) < 2.1)
    {
      err = 0.009;
      applySys[1] = 1;
    }
  else
    {
      err = 0.027;
      applySys[2] = 1;
    }
  
  for(int i=0; i<3; i++){
    string histName = histNumber + "_CMS_scale_m_" +regions[i] +"_2017";
    if(applySys[i] == 0) err = 0;
    weight = event_weight * (1 + err);
    fillHist(histName+"_up", MuIndex, TauIndex, true, weight);
    weight = event_weight * (1 - err);
    fillHist(histName+"_down", MuIndex, TauIndex, true, weight);
  }
    
}
void mutau_analyzer::muonMissID( string histNumber , int eleIndex, int tauIndex, bool isFakeBkg, float event_weight){

  if( !(my_genmatching_l2==2 || my_genmatching_l2==4 || my_genmatching_l1==2 || my_genmatching_l1==4) )
    return;
  
  string region = ""; string taudm = "";
  float err = 0.01;
  double weight = 1.0;
  
  for(int i=0; i<2; i++){
    string histName = histNumber + "_CMS_scale_mfaket_" + tauDM_map[i].first +"_2017";
    if(tauDM_map[i].second != tau_DecayMode->at(tauIndex)) err = 0;
    weight = event_weight * (1 + err);
    fillHist(histName+"_up", MuIndex, TauIndex, true, weight);
    weight = event_weight * (1 - err);
    fillHist(histName+"_down", MuIndex, TauIndex, true, weight);
  }
}

void mutau_analyzer::tauIDunc( string histNumber , int eleIndex, int tauIndex, bool isFakeBkg, float event_weight){
  double sf = tauSFTool.getSFvsPT(my_tauP4.Pt(), my_genmatching_l2);
  double sf_up   = tauSFTool.getSFvsPT(my_tauP4.Pt(), my_genmatching_l2, "Up");
  double sf_down = tauSFTool.getSFvsPT(my_tauP4.Pt(), my_genmatching_l2, "Down");
  //cout<<"tau id sf unc :  "<<sf<<"  "<<sf_up<<"  "<<sf_down<<endl;
  string regions[3] = { "CMS_tauideff_pt30to35_2017", "CMS_tauideff_pt35to40_2017", "CMS_tauideff_ptgt40_2017" };
  bool applySys[3] = {0, 0, 0};
  
  if (my_tauP4.Pt() >= 30 && my_tauP4.Pt()<35)
    applySys[0] = 1;
  else if (my_tauP4.Pt() >= 35 && my_tauP4.Pt()<40)
    applySys[1] = 1;
  else if (my_tauP4.Pt() >= 40 )
    applySys[2] = 1;

  for(int i=0; i<3; i++){
    double weight_up = sf_up/sf;
    double weight_down = sf_down/sf;
    if (applySys[i] == 0) { weight_up = 1; weight_down = 1; }
    fillHist(histNumber+"_"+regions[i]+"_up", MuIndex, TauIndex, true, event_weight*sf_up/sf);
    fillHist(histNumber+"_"+regions[i]+"_down", MuIndex, TauIndex, true, event_weight*sf_down/sf);
  }

}
void mutau_analyzer::dyShape( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  
  fillHist(histNumber+"_CMS_htt_dyShape_up", MuIndex, TauIndex, true, event_weight * 1.1);
  fillHist(histNumber+"_CMS_htt_dyShape_down", MuIndex, TauIndex, true, event_weight * 0.9);
}
void mutau_analyzer::ttbarShape( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  
  double top_pt_weight=1.0;
  if(found_TTbar_sample){
    int t_index = get_t_Cand(); int tbar_index = get_tbar_Cand();
    if( t_index >-1 && tbar_index > -1 ){
      top_pt_weight = sqrt( TMath::Exp( 0.0615 - 0.0005*mcPt->at(t_index)) * TMath::Exp( 0.0615 - 0.0005*mcPt->at(tbar_index)) );
    }
  }
  fillHist(histNumber+"_CMS_htt_ttbarShape_up", MuIndex, TauIndex, true, event_weight * (2*(top_pt_weight -1 ) + 1) );
  fillHist(histNumber+"_CMS_htt_ttbarShape_down", MuIndex, TauIndex, true, event_weight / top_pt_weight );
}
void mutau_analyzer::metSysUnc( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  
  TLorentzVector mymetvector = my_metP4;
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
  float pfmet=my_metP4.Pt(); float pfmetPhi=my_metP4.Phi();
  float pfmetcorr_ex_responseUp=pfmet*cos(pfmetPhi);     float pfmetcorr_ey_responseUp=pfmet*sin(pfmetPhi);
  float pfmetcorr_ex_responseDown=pfmet*cos(pfmetPhi);   float pfmetcorr_ey_responseDown=pfmet*sin(pfmetPhi);
  float pfmetcorr_ex_resolutionUp=pfmet*cos(pfmetPhi);   float pfmetcorr_ey_resolutionUp=pfmet*sin(pfmetPhi);
  float pfmetcorr_ex_resolutionDown=pfmet*cos(pfmetPhi); float pfmetcorr_ey_resolutionDown=pfmet*sin(pfmetPhi);
  float genpX = BosonP4.Px();
  float genpY = BosonP4.Py();
  float vispX = visGenP4.Px();
  float vispY = visGenP4.Py();
  int recoiljets = my_njets;
  int Process= MEtSys::ProcessType::BOSON;
  metSys.ApplyMEtSys(mymetvector.Px(),mymetvector.Py(),
		     genpX,genpY,
		     vispX,vispY,
		     recoiljets,
		     Process,
		     MEtSys::SysType::Response,
		     MEtSys::SysShift::Up,
		     pfmetcorr_ex_responseUp,pfmetcorr_ey_responseUp
		     );

  metSys.ApplyMEtSys(mymetvector.Px(),mymetvector.Py(),
		     genpX,genpY,
		     vispX,vispY,
		     recoiljets,
		     Process,
		     MEtSys::SysType::Response,
		     MEtSys::SysShift::Down,
		     pfmetcorr_ex_responseDown,pfmetcorr_ey_responseDown
		     );

  metSys.ApplyMEtSys(mymetvector.Px(),mymetvector.Py(),
		     genpX,genpY,
		     vispX,vispY,
		     recoiljets,
		     Process,
		     MEtSys::SysType::Resolution,
		     MEtSys::SysShift::Up,
		     pfmetcorr_ex_resolutionUp,pfmetcorr_ey_resolutionUp
		     );
  metSys.ApplyMEtSys(mymetvector.Px(),mymetvector.Py(),
		     genpX,genpY,
		     vispX,vispY,
		     recoiljets,
		     Process,
		     MEtSys::SysType::Resolution,
		     MEtSys::SysShift::Down,
		     pfmetcorr_ex_resolutionDown,pfmetcorr_ey_resolutionDown
		     );
  // cout<<" ********************************************"<<endl;
  // cout<<"nominal    "<<" "<<mymetvector.Px()<<"  "<<mymetvector.Py()<<endl;
  // cout <<"response"<<"\n"<<" "<<pfmetcorr_ex_responseUp<<" "<<pfmetcorr_ex_responseDown<<endl;
  // cout <<"response"<<"\n"<<" "<<pfmetcorr_ey_responseUp<<" "<<pfmetcorr_ey_responseDown<<endl;
  // cout <<"resolution"<<"\n"<<" "<<pfmetcorr_ex_resolutionUp<<" "<<pfmetcorr_ex_resolutionDown<<endl;
  // cout <<"resolution"<<"\n"<<" "<<pfmetcorr_ey_resolutionUp<<" "<<pfmetcorr_ey_resolutionDown<<endl;
  string njetName ="";
  // if (recoiljets ==0) njetName="0";
  // else if (recoiljets ==1) njetName="1";
  // else if (recoiljets >=2) njetName="2";
  TLorentzVector mymet_responseUp, mymet_responseDown, mymet_resolutionUp, mymet_resolutionDown ;
  mymet_responseUp.SetPxPyPzE(pfmetcorr_ex_responseUp,pfmetcorr_ey_responseUp,0,sqrt(pfmetcorr_ex_responseUp*pfmetcorr_ex_responseUp+pfmetcorr_ey_responseUp*pfmetcorr_ey_responseUp));
  my_metP4 = mymet_responseUp;
  fillHist(histNumber+"_CMS_htt_boson_scale_met_"+njetName+"Jet_2017_up", MuIndex, TauIndex, true, event_weight );
  mymet_responseDown.SetPxPyPzE(pfmetcorr_ex_responseDown,pfmetcorr_ey_responseDown,0,sqrt(pfmetcorr_ex_responseDown*pfmetcorr_ex_responseDown+pfmetcorr_ey_responseDown*pfmetcorr_ey_responseDown));
  my_metP4 = mymet_responseDown;
  fillHist(histNumber+"_CMS_htt_boson_scale_met_"+njetName+"Jet_2017_down", MuIndex, TauIndex, true, event_weight );
  mymet_resolutionUp.SetPxPyPzE(pfmetcorr_ex_resolutionUp,pfmetcorr_ey_resolutionUp,0,sqrt(pfmetcorr_ex_resolutionUp*pfmetcorr_ex_resolutionUp+pfmetcorr_ey_resolutionUp*pfmetcorr_ey_resolutionUp));
  my_metP4 = mymet_resolutionUp;
  fillHist(histNumber+"_CMS_htt_boson_reso_met_"+njetName+"Jet_2017_up", MuIndex, TauIndex, true, event_weight );
  mymet_resolutionDown.SetPxPyPzE(pfmetcorr_ex_resolutionDown,pfmetcorr_ey_resolutionDown,0,sqrt(pfmetcorr_ex_resolutionDown*pfmetcorr_ex_resolutionDown+pfmetcorr_ey_resolutionDown*pfmetcorr_ey_resolutionDown));
  my_metP4 = mymet_resolutionDown;
  fillHist(histNumber+"_CMS_htt_boson_reso_met_"+njetName+"Jet_2017_down", MuIndex, TauIndex, true, event_weight );

  /// reset my_metP4
  my_metP4 = mymetvector;
  
}

void mutau_analyzer::metClusteredUnc( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){

  // if (isData || isSignal) return; ////// why
  TLorentzVector mymet_vector = my_metP4;
  float tmp_pfMET = pfMET;
  float tmp_pfMETPhi = pfMETPhi;
  
  pfMET = pfMET_T1UESUp;
  pfMETPhi = pfMETPhi_T1UESUp;
  resetTauTau( muIndex, tauIndex);
  fillHist(histNumber+"_CMS_scale_met_unclustered_2017_up", MuIndex, TauIndex, true, event_weight);

  pfMET = pfMET_T1UESDo;
  pfMETPhi = pfMETPhi_T1UESDo;  
  resetTauTau( muIndex, tauIndex);
  fillHist(histNumber+"_CMS_scale_met_unclustered_2017_down", MuIndex, TauIndex, true, event_weight);
  
  /////////// we need to reset the jet, met variable
  pfMET = tmp_pfMET;
  pfMETPhi = tmp_pfMETPhi;
  resetTauTau( muIndex, tauIndex);

}

void mutau_analyzer::fillHist( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  string hNumber = histNumber;
  
  plotFill("muPt_"+hNumber,  my_muP4.Pt() , 30 , 20.0 , 80.0,  event_weight);
  plotFill("muEta_"+hNumber, my_muP4.Eta(), 48, -2.4, 2.4,  event_weight);
  plotFill("muPhi_"+hNumber, my_muP4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("muDz_"+hNumber,  muDz->at(muIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("muD0_"+hNumber,  muD0->at(muIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("muonID_"+hNumber, muIDbit->at(muIndex)>>1&1, 4, -2, 2,  event_weight); // electronID
  float relMuIso = ( muPFChIso->at(muIndex) + max( muPFNeuIso->at(muIndex) + muPFPhoIso->at(muIndex) - 0.5 *muPFPUIso->at(muIndex) , 0.0 )) / (my_muP4.Pt());
  plotFill("relMuIso_"+hNumber, relMuIso, 15, 0, 0.3,  event_weight);
  plotFill("muCharge_"+hNumber, muCharge->at(muIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  my_tauP4.Pt() , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, my_tauP4.Eta(), 25, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, my_tauP4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = my_muP4.DeltaR(my_tauP4);
  plotFill("deltaR_"+hNumber, deltaR , 40, 0, 6,  event_weight);
  double deltaPhi = DeltaPhi(muPhi->at(muIndex), tau_Phi->at(tauIndex));
  double deltaEta = fabs(muEta->at(muIndex) - tau_Eta->at(tauIndex));
  plotFill("deltaPhi_"+hNumber, deltaPhi , 30, -3.14, 3.14,  event_weight);
  plotFill("deltaEta_"+hNumber, deltaEta ,  25, -2.5, 2.5,  event_weight);
  
  plotFill("nVertex_"+hNumber, nVtx ,  24, 0, 60,  event_weight);
  plotFill("nJet_"+hNumber,  my_njets , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, my_metP4.Pt() , 20, 0, 200,  event_weight);
  plotFill("metLongXaxis_"+hNumber, my_metP4.Pt() , 10, 100, 200,  event_weight);
  plotFill("metPhi_"+hNumber, my_metP4.Phi() , 30, -3.14, 3.14,  event_weight);
  double mT_muMet = TMass_F( my_muP4.Pt(),my_muP4.Phi(), my_metP4.Pt(), my_metP4.Phi() );
  plotFill("mT_muMet_"+hNumber, mT_muMet , 20, 0, 200,event_weight);

  double visMass_mutau = (my_muP4+ my_tauP4).M();
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  double HiggsPt = (my_muP4+my_tauP4+my_metP4).Pt();
  plotFill("higgsPt_"+hNumber,HiggsPt ,  40, 0, 400,  event_weight);

  double tot_tr_mass = (my_muP4 + my_tauP4 + my_metP4 ).M();
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);
  float TrMassBins[15]={0.0, 20, 40, 60, 90, 120, 150, 180, 210, 235, 260, 285, 325, 400, 2000};
  plotFill_customBinning("tot_TMass_full_"+hNumber, tot_tr_mass , 14, TrMassBins,  event_weight);

  int triggerBin1, triggerBin2, triggerBin3, triggerBin4;
  triggerBin1=triggerBin2=triggerBin3=triggerBin4=0;
  if( HLTEleMuX>>21&1 == 1 )  triggerBin3=3;
  if( HLTEleMuX>>60&1 == 1 ) triggerBin2=2;
  if( HLTTau>>0&1 == 1 )     triggerBin1=1;
  if(triggerBin1>0)
    plotFill("trigger_"+hNumber, triggerBin1 , 5, 0, 5,  event_weight);
  if(triggerBin2>0)
    plotFill("trigger_"+hNumber, triggerBin2 , 5, 0, 5,  event_weight);
  if(triggerBin3>0)
    plotFill("trigger_"+hNumber, triggerBin3 , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(is_MC){
    if(my_genmatching_l2==1) genMatchBin=1;
    else if(my_genmatching_l2==2) genMatchBin=2;
    else if(my_genmatching_l2==3) genMatchBin=3;
    else if(my_genmatching_l2==4) genMatchBin=4;
    else if(my_genmatching_l2==5) genMatchBin=5;
    else if(my_genmatching_l2==6) genMatchBin=6;
  }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}


void mutau_analyzer::makeTestPlot( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight){
  // string hNumber = histNumber;
  // std::vector<int> tmpCand; tmpCand.clear();
  // for(int iEle=0;iEle<nEle;iEle++)
  //   {
  //     tmpCand.push_back(iEle);
  //   }
  // plotFill("elePt_"+hNumber,  elePt->at(tmpCand[0]) , 38 , 24 , 100,  event_weight);
  //cout<<"     elePt_"<<hNumber<<" = "<< elePt->at(tmpCand[0])<<endl;
}


TLorentzVector mutau_analyzer::MetRecoilCorrections(int muIndex, int tauIndex, TLorentzVector mymet){
  //// met recoil correction
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
  float pfmet=mymet.Pt(); float pfmetPhi=mymet.Phi();
  float pfmetcorr_ex=pfmet*cos(pfmetPhi); float pfmetcorr_ey=pfmet*sin(pfmetPhi);
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(muIndex, tauIndex);
  
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

void mutau_analyzer::applyTauESCorrections(TLorentzVector tauP4, int tauIndex, TLorentzVector& tauP4Corr)
{
  int tZTTGenMatching = myGenMaching(tauIndex);
  if(is_MC)
  {
    if (tZTTGenMatching>=5 && tau_DecayMode->at(tauIndex)==0) tauP4Corr=tauP4*1.007;
    else if (tZTTGenMatching>=5 && tau_DecayMode->at(tauIndex)==1) tauP4Corr=tauP4*0.998;
    else if (tZTTGenMatching>=5 && tau_DecayMode->at(tauIndex)==10) tauP4Corr=tauP4*1.001;
    if (  (tZTTGenMatching==1 || tZTTGenMatching==3) && tau_DecayMode->at(tauIndex)==0 ) tauP4Corr=tauP4*0.982;
    else if ( (tZTTGenMatching==1 || tZTTGenMatching==3) && tau_DecayMode->at(tauIndex)==1) tauP4Corr=tauP4*1.018;
    if ( (tZTTGenMatching==2 || tZTTGenMatching==4) && tau_DecayMode->at(tauIndex)==0) tauP4Corr=tauP4*0.998;
    else if ( (tZTTGenMatching==2 || tZTTGenMatching==4) && tau_DecayMode->at(tauIndex)==1) tauP4Corr=tauP4*0.992;
    //tauP4Corr = tauP4Corr* get_zptmass_weight(); 
  }
  else
    tauP4Corr = tauP4;

}


int mutau_analyzer::if_DY_Genmatching(int muIndex, int tauIndex){
  // if(!is_MC)
  //   return 1;
  if( found_DYjet_sample==false )
    return 1;
  else if(found_DYjet_sample==true){
    if(  my_genmatching_l2<5 &&  my_genmatching_l1<5 ) // dy -> ll genmatched
      return 2;
    if (  my_genmatching_l2>=5 &&  my_genmatching_l1<5 ) // dy -> ltau genmatched
      return 3;
  }
  return 0;

}

int mutau_analyzer::eventCategory(int muIndex, int tauIndex, double higgsPt){
  int category=0;
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(muIndex, tauIndex);
  int njets = jetCand.size();
  double mjj=0;
  TLorentzVector jet1P4, jet2P4;
  if(njets>=2)
    {
      jet1P4.SetPtEtaPhiE(jetPt->at(jetCand[0]), jetEta->at(jetCand[0]), jetPhi->at(jetCand[0]), jetE->at(jetCand[0]));
      jet2P4.SetPtEtaPhiE(jetPt->at(jetCand[1]), jetEta->at(jetCand[1]), jetPhi->at(jetCand[1]), jetE->at(jetCand[1]));
      mjj = (jet1P4+ jet2P4).M();
    }
  ///////category selection
  if(njets==0)
    {
      if( higgsPt<10)
        return category=1;
      else if (higgsPt>10)
        return category=2;
    }
  if(njets>=2 && mjj > 350)
    {
      if      (higgsPt<200)
        return category=5;
      else if (higgsPt>200)
        return category=6;
    }
  if(njets==1)
    return category=3;
  if(njets>=2 && mjj<350)
    return category=4;
}
double mutau_analyzer::getTauFES(int tauIndex, string shift){
  int binNumber = 0;
  if(tau_DecayMode->at(tauIndex)==0)
    {
      if (abs(my_tauP4.Eta())< 1.5) 
	binNumber = 0;
      else
	binNumber = 2;
    }
  else if(tau_DecayMode->at(tauIndex)==1)
    {
      if (abs(my_tauP4.Eta())< 1.5) 
	binNumber = 1;
      else
	binNumber = 3;
    }
  float yValue = h_taufesSF->GetY()[binNumber];
  if(shift=="")          return yValue;
  else if(shift=="up")   return yValue + h_taufesSF->GetErrorYhigh(binNumber);
  else if(shift=="down") return yValue - h_taufesSF->GetErrorYlow(binNumber);
  return 1.0;
  
}

double mutau_analyzer::get_zptmass_weight(){
  double weight = 1.0;
  int genZCand= -1;
  if(is_MC)
    genZCand = getZCand();
  if(genZCand>-1){
    w->var("m_pt")->setVal(my_muP4.Pt());
    w->var("m_eta")->setVal(my_muP4.Eta());
    w->var("z_gen_pt")->setVal(mcPt->at(genZCand));
    w->var("z_gen_mass")->setVal(mcMass->at(genZCand));
    weight=w->function("zptmass_weight_nom")->getVal();
    
  }
  return weight;

}
double mutau_analyzer::btag_sf_weight(int muIndex , int tauIndex){
  double weight = 1.0;
  vector<int> looseBjet = bJet_loose(muIndex , tauIndex);
  vector<int> mediumBjet = bJet_medium(muIndex , tauIndex);
  if( mediumBjet.size()>0 )
    {
      weight =  btag_csv->EvalSF(0,"comb","central",1
				  ,jetPt->at(mediumBjet[0])
				  ,jetEta->at(mediumBjet[0])
				  );
      return weight;
    }
  if( looseBjet.size()>0  )
    {
      weight =  btag_csv->EvalSF(0,"comb","central",1
				 ,jetPt->at(looseBjet[0])
				 ,jetEta->at(looseBjet[0])
				 );
      return weight;
    }
  return 1.0;
}
