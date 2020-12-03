////tautau_analyzer.C
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
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"

#include "commonFunctions.h"
#include "fractions.C"
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
  
  tautau_analyzer t(argv[1],argv[2], isMC, SampleName);
  t.Loop(maxEvents,reportEvery, SampleName );
  //delete myMap1;
  cout<<" Outpt written to "<<outputfile<<endl;
  sw.Stop();
  sw.Print();
  return 0;
}

void tautau_analyzer::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
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
  if (fChain == 0) return;
  //std::vector<int> eleCand;        eleCand.clear();
  std::vector<int> tauCand;        tauCand.clear();
  std::vector<int> tau2Cand;        tau2Cand.clear();
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
       	      
      //eleCand.clear();
      tauCand.clear();   tau2Cand.clear();
       
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
       
      double pileup_sf = 1.0;
      bool passTauTrigger=false;
      int report_i=0;
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
       

       
      ////Trigger bit selection
      // 5 : HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v
      // 6 : HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v
      // 7 : HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v
      if(  HLTTau>>5&1==1 || HLTTau>>6&1==1 || HLTTau>>7&1==1  )
	passTauTrigger=true;
       
      /////
      if(debug)cout<<"entry # : "<<jentry<<endl;
       
      if(debug)cout<<"reco selections begin"<<endl;
      tauCand.clear(); tau2Cand.clear();
      ////// reco selection begin
      if(debug)cout<<"signal region DY->ll -  isolated begin"<<endl;
      ////// DY Z-> ll signal region -  isolated begin
      bool dy_ll_genmatching=false;
      //dy_ll_genmatching=true;

      if(!is_MC)
	event_weight=1.0;
      else
	event_weight=weight;
      if(metFilters==0 )
	{
	  if(debug)cout<<"metfilters selected"<<endl;
	  if(is_MC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	  nMETFiltersPassed_dyll+=event_weight;
	  //makeTestPlot("a_dyll", 0,0,0,event_weight);
	  if(debug)cout<<"genweight applied"<<endl;
	  if( passTauTrigger )
	    {
	      nSingleTrgPassed_dyll+=event_weight;
	      if(debug)cout<<"trigger selected"<<endl;
	      //makeTestPlot("b_dyll", 0,0,0,event_weight);
	      tauCand = getTauCand(40.0,2.1);  ///// tau selected
	      if( tauCand.size() >0 ) 
		{ 
		  nGoodMuonPassed_dyll+=event_weight;
		  if(debug)cout<<"this worked Line 443"<<endl;
		  //makeTestPlot("c_dyll", 0,0,0,event_weight);
		   
		  //tauCand = getTauCand(40.0,2.1);
		  if( tauCand.size()>=2  ) 
		    {
		      nGoodTauPassed_dyll+=event_weight;
		      if(debug)cout<<"this worked Line 424"<<endl;
		      //makeTestPlot("d_dyll", 0,0,0,event_weight);
		       
		      setMyEleTau(tauCand[0], tauCand[1]); 
		      // from here we can use my_eleP4, my_tauP4, my_metP4, etc
		      
		      if( eVetoZTTp001dxyz(Tau1Index, Tau2Index)
			  && mVetoZTTp001dxyz(Tau1Index, Tau2Index)
			  ) Ztt_selector=true;
		      else Ztt_selector=false;
		      
		      if ( TriggerSelection(my_tau1P4, my_tau2P4) )
			{
			  if(Ztt_selector) 
			    {
			      if ( tau_Charge->at(Tau1Index)*tau_Charge->at(Tau2Index)<0  )
				   
				{
				  nGoodMuTauPassed_dyll+=event_weight;
				  //makeTestPlot("e_dyll", 0,0,0,event_weight);
				  
				  if ( MatchTriggerFilter(Tau1Index, Tau2Index) )
				    {
			       
				      if(debug)cout<<"this worked Line 314, SR opp charge passed"<<endl;
				  
				      applySf=1.0;
				      if(is_MC)
				      	applySf=  getScaleFactors( my_tau1P4.Pt(),
				      				   my_tau2P4.Pt(),
				      				   my_tau1P4.Eta(),				      				   
				      				   my_tau2P4.Eta(),
				      				   tau_DecayMode->at(Tau1Index),
				      				   myGenMaching(Tau1Index),
				      				   false  /// this is set to true for fake bakground
				      				   );
				   
				      // if(debug)cout<<" sf : "<<getScaleFactors( EleIndex[0] , TauIndex[0] , false , is_MC , debug ) <<endl;
				      event_weight = event_weight * applySf;
				      
				      //makeTestPlot("f_dyll", 0,0,0,event_weight);
				      if( thirdLeptonVeto(Tau1Index, Tau2Index)  )
					{
					  nPassedThirdLepVeto_dyll+=event_weight;
					  
					  //if( pass_bjet_veto )
					  {
					  
					    nPassedBjetVeto_dyll+=event_weight;
					    //makeTestPlot("h_dyll", 0,0,0,event_weight);
					    double deltaR =  my_tau1P4.DeltaR(my_tau2P4);
					    if(deltaR > 0.5 )
					      {
						nDeltaRPassed_dyll+=event_weight;
						if(is_MC==false)event_weight=1.0;
						//makeTestPlot("i_dyll", 0,0,0,event_weight);
						if(debug)cout<<"this worked Line 374"<<endl;
						fillHist("5_dyll",  Tau1Index, Tau2Index, false, event_weight);
					       
						double mT_tauMet = TMass_F(my_tau2P4.Pt(),my_tau2P4.Phi(), my_metP4.Pt(), my_metP4.Phi() );
				
						if( mT_tauMet < 50 )
						  {
						    fillHist("6_dyll", Tau1Index, Tau2Index, false, event_weight);
						    //makeTestPlot("j_dyll", 0,0,0,event_weight);
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
      
      if(debug)cout<<"signal region -  isolated begin L523"<<endl;       
      
      Ztt_selector=false;
      ////// signal region -  isolated begin
      if(is_MC)
	event_weight=weight;
      else
	event_weight=1.0;
      tauCand.clear();   tau2Cand.clear();
	   
      if(metFilters==0)
	{
	  if(debug)cout<<"metfilters selected"<<endl;
	  if (is_MC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	  nMETFiltersPassed+=event_weight;
	  //makeTestPlot("a", 0,0,0,event_weight);
	  if(debug)cout<<"genweight applied"<<endl;
	  if( passTauTrigger  )
	    {
	      nSingleTrgPassed+=event_weight;
	      if(debug)cout<<"trigger selected"<<endl;
	      //makeTestPlot("b", 0,0,0,event_weight);
	      tauCand = getTauCand(40.0,2.1);  ///// tau selected
	      if( tauCand.size() >0 ) 
		{ 
		  nGoodMuonPassed+=event_weight;
		  if(debug)cout<<"this worked Line 526"<<endl;
		       
		  //makeTestPlot("c", 0,0,0,event_weight);
		  if( tauCand.size() >=2 )
		    {
		      nGoodTauPassed+=event_weight;
		      setMyEleTau(tauCand[0], tauCand[1]);
		      //makeTestPlot("d", 0,0,0,event_weight);
			   
		      if( eVetoZTTp001dxyz(Tau1Index, Tau2Index)
                          && mVetoZTTp001dxyz(Tau1Index, Tau2Index)
                          ) Ztt_selector=true;
                      else Ztt_selector=false;			   
			   
		      if ( TriggerSelection(my_tau1P4, my_tau2P4) )
                        {
			  if(Ztt_selector) 
			    {
			   		   
			      if ( tau_Charge->at(Tau1Index)*tau_Charge->at(Tau2Index)<0  )
				
				{
				  nGoodMuTauPassed+=event_weight;
			       
				  if(debug)cout<<"this worked Line 538"<<endl;
			       
				  //makeTestPlot("e", 0,0,0,event_weight);
				  if ( MatchTriggerFilter(Tau1Index, Tau2Index) )
				    {
				      if(debug)cout<<"this worked Line 534"<<endl;
				      applySf=1.0;
				      if(is_MC)
					applySf=  getScaleFactors( my_tau1P4.Pt(),
                                                                   my_tau2P4.Pt(),
                                                                   my_tau1P4.Eta(),
                                                                   my_tau2P4.Eta(),
                                                                   tau_DecayMode->at(Tau1Index),
                                                                   myGenMaching(Tau1Index),
                                                                   false  /// this is set to true for fake bakground
                                                                   );
				      // if(debug)cout<<" sf : "<<getScaleFactors( EleIndex[0] , TauIndex[0] , false , is_MC , debug ) <<endl;
				      event_weight = event_weight * applySf;
				      //makeTestPlot("f", 0,0,0,event_weight);
				      if( thirdLeptonVeto(Tau1Index, Tau2Index)  )
					{
					  nPassedThirdLepVeto+=event_weight;
					  //makeTestPlot("g", 0,0,0,event_weight);
					  //if( pass_bjet_veto )
					  {
					    nPassedBjetVeto+=event_weight;
					    //makeTestPlot("h", 0,0,0,event_weight);
					    //if(tau_DecayMode->at(TauIndex)==5 || tau_DecayMode->at(TauIndex)==6) continue;
					    
					    double deltaR = my_tau1P4.DeltaR(my_tau2P4);
					    if(deltaR > 0.5 )
					      {
						nDeltaRPassed+=event_weight;
						if(is_MC==false)event_weight=1.0;
						if(debug)cout<<"this worked Line 558"<<endl;
						fillHist("5",  Tau1Index, Tau2Index, false, event_weight);
						//makeTestPlot("i", 0,0,0,event_weight);
						double mT_tauMet = TMass_F(my_tau2P4.Pt(),my_tau2P4.Phi(), my_metP4.Pt(), my_metP4.Phi() );
						  
						if( mT_tauMet < 50 )
						  {
						    fillHist("6", Tau1Index, Tau2Index, false, event_weight);
						    //makeTestPlot("j", 0,0,0,event_weight);
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
      if(is_MC)
	event_weight=weight;
      else
	event_weight=1.0;
      tauCand.clear();   tau2Cand.clear();
      if(metFilters==0)
	{
	  if (is_MC)fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	  nMETFiltersPassed_fr+=event_weight;
	  //makeTestPlot("a_fr", 0,0,0,event_weight);
	  if( passTauTrigger  )
	    {
	      nSingleTrgPassed_fr+=event_weight;
	      if(debug)cout<<"trigger selected line 636"<<endl;
	      //makeTestPlot("b_fr", 0,0,0,event_weight);
	      bool tautauSelected=false;
	      int leadingtau = -1; int subleadingtau = -1;
	      tauCand = getTauCand_noID(40.0,2.1);  ///// taus selected
	      if( tauCand.size()>0 )
		{
		  nGoodMuonPassed_fr+=event_weight;
		  //makeTestPlot("c_fr", 0,0,0,event_weight);
		  if(debug)cout<<"this worked Line "<<__LINE__<<endl;
		  assignTauIndices(tauCand , leadingtau, subleadingtau);  // select iso and anti-iso taus -> subleading & leading taus
		  if( tauCand.size()>=2 && leadingtau>-1 && subleadingtau>-1 ) 
		    { 
		      // cout<<"getTauCand_noID size ="<<tauCand.size()
		      // 	  <<"  leadingtau :"<< tau_Pt->at(leadingtau)
		      // 	  <<"  subleadingtau : " << tau_Pt->at(subleadingtau)
		      // 	  << " leadingtau :"<<tau_byVVVLooseDeepTau2017v2p1VSjet->at(leadingtau)
		      // 	  << tau_byMediumDeepTau2017v2p1VSjet->at(leadingtau)
		      // 	  << " subleading tau : " << tau_byMediumDeepTau2017v2p1VSjet->at(subleadingtau)
		      // 	  <<endl;
		      
		      //nGoodTauPassed_fr+=event_weight;
		      //makeTestPlot("d_fr", 0,0,0,event_weight);
		      
		      setMyEleTau(leadingtau, subleadingtau); 
		      double mt=TMass_F(my_tau2P4.Pt(),my_tau2P4.Phi()
					,my_metP4.Pt(), my_metP4.Phi());
		      double mvis=(my_tau1P4+my_tau2P4).M();
		      double higgsPt = (my_tau1P4 + my_tau2P4 + my_metP4).Pt();
		      double frac_tt=0.01; double frac_qcd=0.24; double frac_w=0.75; 
		      int category=eventCategory(Tau1Index , Tau2Index, higgsPt) ;
		      getFractions(category, mvis, frac_qcd, frac_w, frac_tt); /// this assigns right values for qcd, w and tt fractions
		      float my_fakefactor = get_ff( my_tau1P4.Pt(), my_tau2P4.Pt(), mt, mvis, my_njets
						    , frac_tt, frac_qcd, frac_w
						    , ff_qcd_0jet, ff_qcd_1jet,ff_qcd_2jet
						    , subleadingTauPt_correction
						    , osssclosure_qcd);
		      //printTabSeparated("my_fakefactor",my_fakefactor,
		      //			"fractions",frac_qcd, frac_w, frac_tt
		      //			);
		      event_weight = event_weight*my_fakefactor;
		      if( eVetoZTTp001dxyz(Tau1Index, Tau2Index)
			  && mVetoZTTp001dxyz(Tau1Index, Tau2Index)
			  ) Ztt_selector=true;
		      else Ztt_selector=false;
		      nGoodTauPassed_fr+=event_weight;
		      if ( TriggerSelection(my_tau1P4, my_tau2P4) )
			{
			  if(Ztt_selector)
			    {
			      if ( tau_Charge->at(Tau1Index)*tau_Charge->at(Tau2Index)<0 )
				{				  
				  nGoodMuTauPassed_fr+=event_weight;
				  //makeTestPlot("e_fr", 0,0,0,event_weight);
				  if ( MatchTriggerFilter(Tau1Index, Tau2Index) )
				    {
				      
				      applySf=1.0;
				      if(is_MC)
                                        applySf=  getScaleFactors( my_tau1P4.Pt(),
                                                                   my_tau2P4.Pt(),
                                                                   my_tau1P4.Eta(),
                                                                   my_tau2P4.Eta(),
                                                                   tau_DecayMode->at(Tau1Index),
                                                                   myGenMaching(Tau1Index),
                                                                   false  /// this is set to true for fake bakground
                                                                   );
				  
				      event_weight = event_weight * applySf;
				      
				      //event_weight = event_weight* getFR(TauIndex);

				      //makeTestPlot("f_fr", 0,0,0,event_weight);
				      if( thirdLeptonVeto(Tau1Index, Tau2Index) )
					{
					  nPassedThirdLepVeto_fr+=event_weight;
					  //makeTestPlot("g_fr", 0,0,0,event_weight);
					  //if( pass_bjet_veto )
					  {
					    nPassedBjetVeto_fr+=event_weight;
					    //makeTestPlot("h_fr", 0,0,0,event_weight);
					    //if(tau_DecayMode->at(TauIndex)==5 || tau_DecayMode->at(TauIndex)==6) continue;
					   
					    double deltaR = my_tau1P4.DeltaR(my_tau2P4);
					    if(deltaR > 0.5 )
					      {
						nDeltaRPassed_fr+=event_weight;
						//makeTestPlot("i_fr", 0,0,0,event_weight);
						if(debug)cout<<"this worked Line 442"<<endl;
						fillHist("5_fr", Tau1Index, Tau2Index, true, event_weight);
						//double mT_tauMet = TMass_F(my_tau2P4.Pt(),my_tau2P4.Phi(), my_metP4.Pt(), my_metP4.Phi() );
						  
						if( mt < 50 )
						    
						  {
						    fillHist("6_fr", Tau1Index, Tau2Index, true, event_weight);
						    //makeTestPlot("j_fr", 0,0,0,event_weight);
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
   
  ///

  // fileName->cd();
  // map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
  // map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
  // for (; iMap1 != jMap1; ++iMap1)
  //   nplot1(iMap1->first)->Write();
}

void tautau_analyzer::BookHistos(const char* file1, const char* file2)
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


void tautau_analyzer::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{
  
  //  h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);

}




//---------------------------------------------------                                                                                                                                
// get a electron candiate based on pt eta and isolation                                                                                                                               
//----------------------------------------------------                                                                                                                               

std::vector<int> tautau_analyzer::getEleCand(double elePtCut, double eleEtaCut){
  std::vector<int> tmpCand;
  tmpCand.clear();
  TLorentzVector dau1;
  //Loop over electrons                                                                     
  for(int iEle=0;iEle<nEle;iEle++)
    {
      dau1.SetPtEtaPhiE(elePt->at(iEle), eleEta->at(iEle), elePhi->at(iEle), eleE->at(iEle));
      //dau1 = dau1*(eleCalibE->at(iEle)/dau1.E());
      applyEleESCorrections(dau1, iEle, dau1);

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
      if( kinematic && electronId && relative_iso ){
	tmpCand.push_back(iEle);
      }	
    }                           
  return tmpCand;
  
}

std::vector<int> tautau_analyzer::getTauCand_noID(double tauPtCut, double tauEtaCut){
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
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool mutau_separation=false;
      bool newDecayModeFinding=false;
      bool tau_reject=false;
      bool trigger = false;
      if( dau2.Pt() > tauPtCut 
	  && fabs( dau2.Eta() )< tauEtaCut 
	  && tau_LeadChargedHadron_dz->at(iTau) < 0.2
	  && fabs(tau_Charge->at(iTau))==1
	  )kinematic = true;
      if( tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 ) tauIsolation=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
      
      if( kinematic==true    
	  && decayModeCut==true   
	  && tau_reject==true   
	  && newDecayModeFinding==true
	  && tauIsolation==true
	  )
	{
	  tmpCand.push_back(iTau);
    	}                                     
    }                                                                                       
  return tmpCand;
  
}
void tautau_analyzer::assignTauIndices(vector<int> tauIndices , int& tau1 , int& tau2 ){
  if( tauIndices.size()<2 )
    {
      tau1=-1; tau2=-1; return;
    }
  else if( tauIndices.size()==2 )
    {      
      if(tau_byVVVLooseDeepTau2017v2p1VSjet->at(tauIndices[0])==1
	 && !(tau_byMediumDeepTau2017v2p1VSjet->at(tauIndices[0])==1)
	 && tau_byMediumDeepTau2017v2p1VSjet->at(tauIndices[1])==1
	 )
 	{ tau1 = tauIndices[0]; tau2 = tauIndices[1]; return; }
    }
  else  if( tauIndices.size()>2 )
    {
      for(int ai=0; ai<tauIndices.size(); ai++)
	{
	  for(int i=0; i<tauIndices.size(); i++)
	    {
	      if( tau_byVVVLooseDeepTau2017v2p1VSjet->at(tauIndices[ai])==1 
		  && !(tau_byMediumDeepTau2017v2p1VSjet->at(tauIndices[ai])==1)
		  && tau_byMediumDeepTau2017v2p1VSjet->at(tauIndices[i])==1
		  && tau_Pt->at(tauIndices[ai]) > tau_Pt->at(tauIndices[i])
		  )
		{
		  tau1 = tauIndices[ai]; tau2 = tauIndices[i];
		  return;
		}
	    }
	}
    }
  return;
}

std::vector<int> tautau_analyzer::getTauCand(double tauPtCut, double tauEtaCut){
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
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool mutau_separation=false;
      bool newDecayModeFinding=false;
      bool tau_reject=false;
      bool trigger = false;
      if( dau2.Pt() > tauPtCut 
	  && fabs( dau2.Eta() )< tauEtaCut 
	  && tau_LeadChargedHadron_dz->at(iTau) < 0.2
	  && fabs(tau_Charge->at(iTau))==1
	  )kinematic = true;
      if( tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1 ) tauIsolation=true; 
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;

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
std::vector<int> tautau_analyzer::getAISRTauCand(double tauPtCut, double tauEtaCut){
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
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;
      bool mutau_separation=false;
      bool newDecayModeFinding=false;
      bool tau_reject=false;
      bool trigger = false;
      if( dau2.Pt() > tauPtCut
          && fabs( dau2.Eta() )< tauEtaCut
          && tau_LeadChargedHadron_dz->at(iTau) < 0.2
          && fabs(tau_Charge->at(iTau))==1
          )kinematic = true;
      if(  tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 && !(tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1) ) tauIsolation=true;
      if( tau_DecayMode->at(iTau)==0 || tau_DecayMode->at(iTau)==1 || tau_DecayMode->at(iTau)==10 || tau_DecayMode->at(iTau)==11 ) decayModeCut=true;
      if( tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau)==1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau)==1)tau_reject=true;
      if( tau_IDbits->at(iTau)>>1&1==1 ) newDecayModeFinding=true;
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
int tautau_analyzer::getZCand()
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
int tautau_analyzer::get_t_Cand()
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
int tautau_analyzer::get_tbar_Cand()
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
std::vector<int> tautau_analyzer::getJetCand(int tau1Index, int tau2Index){
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
      
      double lepton1Phi=tau_Phi->at(tau1Index);
      double lepton1Eta= tau_Eta->at(tau1Index);
      double lepton2Phi=0;double lepton2Eta=0;
      lepton2Phi= tau_Phi->at(tau2Index); lepton2Eta=tau_Eta->at(tau2Index);
      double dr_jetEle=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton1Phi, lepton1Eta );
      double dr_jetTau=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton2Phi, lepton2Eta);
      if( dr_jetEle>0.5 && dr_jetTau>0.5 )
        drPassed=true;

      if(kinematic30 && !foundNoisyJets && passLoosePUID && drPassed==true)
        tmpCand.push_back(iJet);
    }
  return tmpCand;
}

int tautau_analyzer::thirdLeptonVeto(){
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
bool tautau_analyzer::thirdLeptonVeto(int tau1Index, int tau2Index){
  std::vector<int> tmpCand; 
  tmpCand.clear();
  int thirdLepIndex = -1;
  std::vector<int> tmpEleCand;
  tmpEleCand.clear();
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
	tmpEleCand.push_back(iEle);
      }                                                                                      
    }                                                          
  double deltaRm1=0; double deltaRm2=0;
  if(tmpEleCand.size() > 0 )
    { 
      deltaRm1 = delta_R(tau_Phi->at(tau1Index),tau_Eta->at(tau1Index), elePhi->at(tmpEleCand[0]),  eleEta->at(tmpEleCand[0]));
      deltaRm2 = delta_R(tau_Phi->at(tau1Index),tau_Eta->at(tau1Index), elePhi->at(tmpEleCand[0]),  eleEta->at(tmpEleCand[0]));
      if(deltaRm1>0.5 || deltaRm2>0.5 ){
	return false;
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
      if(muonId==true && kinematic==true && relative_iso==true){
	tmpCand.push_back(iMu);
      }                   
    } 
  deltaRm1=0;  deltaRm2=0;
  if(tmpCand.size() > 0 )
    { 
      deltaRm1 = delta_R(tau_Phi->at(tau1Index),tau_Eta->at(tau1Index), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      deltaRm2 = delta_R(tau_Phi->at(tau2Index),tau_Eta->at(tau2Index), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      if(deltaRm1>0.5 || deltaRm2>0.5 ){
	return false;
      }
    }
  
  return true;
  
}
bool tautau_analyzer::thirdLeptonVeto(int eleIndex, int tauIndex, int ele2Index){
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
float tautau_analyzer::pTvecsum_F(TLorentzVector a, TLorentzVector b, TLorentzVector c) {
  float pt_vecSum = (a + b+ c).Pt();
  return pt_vecSum;
}

vector<int> tautau_analyzer::bJet_medium(int tau1Index, int tau2Index)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  double lepton1Phi=tau_Phi->at(tau1Index);
  double lepton1Eta= tau_Eta->at(tau1Index);
  double lepton2Phi= tau_Phi->at(tau2Index); double lepton2Eta=tau_Eta->at(tau2Index);   
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

vector<int> tautau_analyzer::bJet_loose(int tau1Index, int tau2Index)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  double lepton1Phi=tau_Phi->at(tau1Index);
  double lepton1Eta= tau_Eta->at(tau1Index);
  double lepton2Phi= tau_Phi->at(tau2Index); double lepton2Eta=tau_Eta->at(tau2Index);   
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
bool tautau_analyzer::passDiElectronVeto(int eleIndex)
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
bool tautau_analyzer::eVetoZTTp001dxyz(int tau1Index, int tau2Index){
  std::vector<int> tmpCand;  tmpCand.clear();
  //std::vector<int> output;  output.clear();
  bool awayFromEverything = true;   int tmpEleIndex=-1;
  //Loop over electrons      
  for(int iEle=0;iEle<nEle;iEle++)
    {
      //if(iEle==eleIndex)continue;
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
      double deltaR_et = delta_R(tau_Phi->at(tau1Index), tau_Eta->at(tau1Index), elePhi->at(tmpCand[0]), eleEta->at(tmpCand[0]));
      double deltaR_ee = delta_R(tau_Phi->at(tau2Index), tau_Eta->at(tau2Index), elePhi->at(tmpCand[0]), eleEta->at(tmpCand[0]));
      if(! (deltaR_et>0.0001 && deltaR_ee>0.0001))
	return false;
    }
  else 
    return true;
    
}
bool tautau_analyzer::mVetoZTTp001dxyz(int tau1Index, int tau2Index){
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
      deltaRm1 = delta_R(tau_Phi->at(tau1Index), tau_Eta->at(tau1Index), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      deltaRm2 = delta_R(tau_Phi->at(tau2Index), tau_Eta->at(tau2Index), muPhi->at(tmpCand[0]),  muEta->at(tmpCand[0]));
      if(! (deltaRm1>0.0001 && deltaRm2>0.0001) )
	return false;
    }
  else
    return true;
  
}
int tautau_analyzer::myGenMaching(int tauIndex)
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
int tautau_analyzer::myGenMaching1(int eleIndex)
{
  if(is_MC==false)
    return 0;
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
std::vector<int> tautau_analyzer::getGenMu(){
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
float tautau_analyzer::exponential(float x,float a,float b,float c) {
  return a * TMath::Exp(-b * x) + c;
}
double tautau_analyzer::exponential( double pT) {
  return TMath::Exp(0.088 - 0.00087*pT + 0.00000092*pow(pT,2) ) ; 
}
double tautau_analyzer::getScaleFactors(  double tau1pt, double tau2pt, double tau1eta, double tau2eta, int taudm, int tau1GenMatch, bool isFakebkg)
{
  bool debug=false;
  double rv_sf=1.0;
  double sf_tau1Trg = 1.0;  double sf_tau2Trg = 1.0;
  double sf_htt_tau1=1.0;   double sf_htt_tau2=1.0;
  double sf_tau1idSF_m = 1.0; double sf_tau2idSF_m = 1.0;
  double sf_tau1idSF_vvvl = 1.0; double sf_tau2idSF_vvvl = 1.0;
  double sf_fakeEle = 1.0; double sf_fakeMu = 1.0;
  double sf_fakeEleES = 1.0; double sf_fakeMuES = 1.0;
  double sf_taufesSF = 1.0;
  
  
  if(  my_genmatching_l1>=5 )
    {
      sf_tau1idSF_m = get_BinContent( h_tauidSF_m, tau1dm) ;
      sf_tau1idSF_vvvl = get_BinContent( h_tauidSF_vvvl, tau1dm);
    }
  if(  my_genmatching_l2>=5 )
    {
      sf_tau2idSF_m = get_BinContent( h_tauidSF_m, tau2dm) ;
      sf_tau2idSF_vvvl = get_BinContent( h_tauidSF_vvvl, tau2dm);
    }
  
  double higgsPt = pTvecsum_F(my_tau1P4, my_tau2P4, my_metP4);
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
      
  
  double tau1PtCheck=tau1pt;
  if(tau1pt > 450 ) tau1PtCheck = 450;
  else if ( tau1pt < 20 )  tau1PtCheck = 20;
  double tau2PtCheck=tau2pt;
  if(tau2pt > 450 ) tau2PtCheck = 450;
  else if ( tau2pt < 20 )  tau2PtCheck = 20;
  //get_BinContent(TH1* histo, double xValue)
  if(tau1dm==0 )  sf_tau1Trg= get_BinContent(h_tauTrgSF_dm0 , tau1PtCheck);
  if(tau1dm==1 )  sf_tau1Trg= get_BinContent(h_tauTrgSF_dm1 , tau1PtCheck);
  if(tau1dm==10)  sf_tau1Trg= get_BinContent(h_tauTrgSF_dm10, tau1PtCheck);
  if(tau1dm==11)  sf_tau1Trg= get_BinContent(h_tauTrgSF_dm11 , tau1PtCheck);

  if(tau2dm==0 )  sf_tau2Trg= get_BinContent(h_tauTrgSF_dm0 , tau2PtCheck);
  if(tau2dm==1 )  sf_tau2Trg= get_BinContent(h_tauTrgSF_dm1 , tau2PtCheck);
  if(tau2dm==10)  sf_tau2Trg= get_BinContent(h_tauTrgSF_dm10, tau2PtCheck);
  if(tau2dm==11)  sf_tau2Trg= get_BinContent(h_tauTrgSF_dm11 , tau2PtCheck);
  ///// btag efficiency
  //double weight_btagSF = btag_sf;
  //cout<<"btag_sf : "<<weight_btagSF<<endl;

 
  //// tau trigger, id scale factors from RooWorkspace
  w->var("t_pt")->setVal(my_tau1P4.Pt());
  w->var("t_eta")->setVal(my_tau1P4.Eta());
  w->var("t_phi")->setVal(my_tau1P4.Phi());
  w->var("t_dm")->setVal(tau1dm);
  double t1_deeptauid = w->function("t_deeptauid_dm_medium")->getVal();
  double t1_trg_pog = w->function("t_trg_pog_deeptau_medium_mutau_ratio")->getVal();
  sf_htt_tau1 = t1_deeptauid*t1_trg_pog;
  
  w->var("t_pt")->setVal(my_tau2P4.Pt());
  w->var("t_eta")->setVal(my_tau2P4.Eta());
  w->var("t_phi")->setVal(my_tau2P4.Phi());
  w->var("t_dm")->setVal(tau2dm);
  double t2_deeptauid = w->function("t_deeptauid_dm_medium")->getVal();
  double t2_trg_pog = w->function("t_trg_pog_deeptau_medium_mutau_ratio")->getVal();
  sf_htt_tau2 = t2_deeptauid*t2_trg_pog;
  

  rv_sf = sf_htt_tau1 * sf_htt_tau2 * sf_tau1Trg * sf_tau2Trg * sf_tau1idSF_m * sf_tau2idSF_m * higgPt_weight ;
  
  if(isFakebkg)
    rv_sf=rv_sf * sf_tau1idSF_vvvl;
  // printTabSeparated("rv_sf :",rv_sf);
  // printTabSeparated("sf_htt_tau : ",sf_htt_tau1 , sf_htt_tau2 );
  // printTabSeparated("higgPt_weight ",higgPt_weight);

  if(rv_sf>0)
    return rv_sf;
  else
    return 1.0;

}
bool tautau_analyzer::TriggerSelection(TLorentzVector tau1P4, TLorentzVector tau2P4){

  if(  HLTTau>>5&1==1 || HLTTau>>6&1==1 || HLTTau>>7&1==1  )
    {
      if(tau1P4.Pt()>40.0 && tau2P4.Pt()>40.0)
        {
          return true;
        }
    }
  return false;
}
bool tautau_analyzer::MatchTriggerFilter(int eleIndex, int tauIndex)
{
  
  return true;
}

double  tautau_analyzer::getFR(int tauIndex){
  double frWeight=1.0;
  return frWeight;
}

void tautau_analyzer::fillHist( string histNumber , int tau1Index, int tau2Index, bool isFakeBkg, float event_weight){
  string hNumber = histNumber;
  /// leading tau 
  plotFill("tau1Pt_"+hNumber,  my_tau1P4.Pt() , 30 , 40 , 180,  event_weight);
  plotFill("tau1Eta_"+hNumber, my_tau1P4.Eta(), 22, -2.1, 2.1,  event_weight);
  plotFill("tau1Phi_"+hNumber, my_tau1P4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("tau1Charge_"+hNumber, tau_Charge->at(tau1Index), 8, -2, 2 ,  event_weight);
  plotFill("tau1Iso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tau1Index), 4, -2, 2,  event_weight);
  plotFill("tau1DecayMode_"+hNumber, tau_DecayMode->at(tau1Index) , 12, 0, 12,  event_weight);
  /// subleading tau
  plotFill("tau2Pt_"+hNumber,  my_tau2P4.Pt() , 30 , 40 , 120,  event_weight);
  plotFill("tau2Eta_"+hNumber, my_tau2P4.Eta(), 22, -2.1, 2.1,  event_weight);
  plotFill("tau2Phi_"+hNumber, my_tau2P4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("tau2Iso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tau2Index), 4, -2, 2,  event_weight);
  plotFill("tau2DecayMode_"+hNumber, tau_DecayMode->at(tau2Index) , 12, 0, 12,  event_weight);
  plotFill("tau2Charge_"+hNumber, tau_Charge->at(tau2Index), 8, -2, 2 ,  event_weight);


  double deltaR = my_tau1P4.DeltaR(my_tau2P4);
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);
  double deltaPhi = DeltaPhi(my_tau1P4.Phi(), my_tau2P4.Phi());
  double deltaEta = fabs(my_tau1P4.Eta() - my_tau2P4.Eta());
  plotFill("deltaPhi_"+hNumber, deltaPhi , 30, 0, 3.14,  event_weight);
  plotFill("deltaEta_"+hNumber, deltaEta ,  25, 0, 2.5,  event_weight);

  plotFill("nJet_"+hNumber,  my_njets , 8, 0, 8,  event_weight);
  plotFill("met_"+hNumber, my_metP4.Pt() , 30, 0, 150,  event_weight);
  plotFill("metLongXaxis_"+hNumber, my_metP4.Pt() , 40, 0, 200,  event_weight);
  plotFill("metPhi_"+hNumber, my_metP4.Phi() , 30, -3.14, 3.14,  event_weight);
  double mT_tauMet = TMass_F( my_tau2P4.Pt(),my_tau2P4.Phi(), my_metP4.Pt(), my_metP4.Phi() );
  plotFill("mT_tauMet_"+hNumber, mT_tauMet , 30, 0, 150,  event_weight);

  double visMass_mutau = (my_tau1P4+ my_tau2P4).M();
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 0, 300,  event_weight);
  
  double HiggsPt = (my_tau1P4+my_tau2P4+my_metP4).Pt();
  plotFill("higgsPt_"+hNumber,HiggsPt , 30, 0, 230,  event_weight);

  double tot_tr_mass = (my_tau1P4 + my_tau2P4 + my_metP4 ).M();
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  int triggerBin1, triggerBin2, triggerBin3, triggerBin4;
  triggerBin1=triggerBin2=triggerBin3=triggerBin4=0;
  if( HLTTau>>7&1 == 1 )     triggerBin3=3;
  if( HLTTau>>6&1 == 1 )     triggerBin2=2;
  if( HLTTau>>5&1 == 1 )     triggerBin1=1;
  if(triggerBin1>0)
    plotFill("trigger_"+hNumber, triggerBin1 , 5, 0, 5,  event_weight);
  if(triggerBin2>0)
    plotFill("trigger_"+hNumber, triggerBin2 , 5, 0, 5,  event_weight);
  if(triggerBin3>0)
    plotFill("trigger_"+hNumber, triggerBin3 , 5, 0, 5,  event_weight);
  if(triggerBin4>0)
    plotFill("trigger_"+hNumber, triggerBin4 , 5, 0, 5,  event_weight);
  
  
}


float tautau_analyzer::EletriggerSF(float pt, float eta){
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
void tautau_analyzer::makeTestPlot( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight){
  string hNumber = histNumber;
  std::vector<int> tmpCand; tmpCand.clear();
  for(int iTau=0;iTau<nTau;iTau++)
    {
      tmpCand.push_back(iTau);
    }
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tmpCand[0]) , 38 , 24 , 100,  event_weight);
  //cout<<"     elePt_"<<hNumber<<" = "<< elePt->at(tmpCand[0])<<endl;
}


TLorentzVector tautau_analyzer::MetRecoilCorrections(int eleIndex, int tauIndex, TLorentzVector mymet){
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
  jetCand=getJetCand(eleIndex, tauIndex);
  
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

void tautau_analyzer::applyTauESCorrections(TLorentzVector tauP4, int tauIndex, TLorentzVector& tauP4Corr)
{
  
  if(is_MC)
  {
    if (myGenMaching(tauIndex)>=5 && tau_DecayMode->at(tauIndex)==0) tauP4Corr=tauP4*1.007;
    else if (myGenMaching(tauIndex)>=5 && tau_DecayMode->at(tauIndex)==1) tauP4Corr=tauP4*0.998;
    else if (myGenMaching(tauIndex)>=5 && tau_DecayMode->at(tauIndex)==10) tauP4Corr=tauP4*1.001;
    if (  (myGenMaching(tauIndex)==1 || myGenMaching(tauIndex)==3) && tau_DecayMode->at(tauIndex)==0 ) tauP4Corr=tauP4*1.003;
    else if ( (myGenMaching(tauIndex)==1 || myGenMaching(tauIndex)==3) && tau_DecayMode->at(tauIndex)==1) tauP4Corr=tauP4*1.036;
    
    //tauP4Corr = tauP4Corr* get_zptmass_weight(); 
  }
  else
    tauP4Corr = tauP4;
}

void tautau_analyzer::applyEleESCorrections(TLorentzVector eleP4, int eleIndex, TLorentzVector& eleP4Corr)
{
  eleP4Corr = eleP4*(eleCalibE->at(eleIndex)/eleP4.E());
  // if(is_MC)
  //   eleP4Corr = eleP4Corr * get_zptmass_weight();
  
}

int tautau_analyzer::if_DY_Genmatching(int eleIndex, int tauIndex){
  // if(!is_MC)
  //   return 1;
  if( found_DYjet_sample==false )
    return 1;
  else if(found_DYjet_sample==true){
    if(  myGenMaching(tauIndex)<5 &&  myGenMaching1(eleIndex)<5 ) // dy -> ll genmatched
      return 2;
    if (  myGenMaching(tauIndex)>=5 &&  myGenMaching1(eleIndex)<5 ) // dy -> ltau genmatched
      return 3;
  }
  return 0;

}

int tautau_analyzer::eventCategory(int tau1Index, int tau2Index, double higgsPt){
  int category=0;
  std::vector<int> jetCand;       jetCand.clear();
  jetCand = getJetCand(tau1Index, tau2Index);
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
double tautau_analyzer::getTauFES(int tauIndex){
  double sf_taufesSF=1.0;
  if(tau_DecayMode->at(tauIndex)==0 && abs(tau_Eta->at(tauIndex))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(1);
  if(tau_DecayMode->at(tauIndex)==0 && abs(tau_Eta->at(tauIndex))>1.4 )  sf_taufesSF = h_taufesSF->Eval(3);
  if(tau_DecayMode->at(tauIndex)==1 && abs(tau_Eta->at(tauIndex))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(5);
  if(tau_DecayMode->at(tauIndex)==1 && abs(tau_Eta->at(tauIndex))>1.4 )  sf_taufesSF = h_taufesSF->Eval(7);
  return sf_taufesSF;

}

double tautau_analyzer::get_zptmass_weight(){
  double weight = 1.0;
  int genZCand= -1;
  if(is_MC)
    genZCand = getZCand();
  if(genZCand>-1){
    w->var("z_gen_pt")->setVal(mcPt->at(genZCand));
    w->var("z_gen_mass")->setVal(mcMass->at(genZCand));
    weight=w->function("zptmass_weight_nom")->getVal();
    
  }
  return weight;

}
double tautau_analyzer::btag_sf_weight(int eleIndex , int tauIndex){
  double weight = 1.0;
  vector<int> looseBjet = bJet_loose(eleIndex , tauIndex);
  vector<int> mediumBjet = bJet_medium(eleIndex , tauIndex);
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
