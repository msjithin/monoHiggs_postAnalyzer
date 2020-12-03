#define smhtt_2017_cxx
#include "smhtt_2017.h"
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
#include "commonFunctions.h"
#include "ApplyFF.h"

using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{
  TStopwatch sw;
  sw.Start();
  
  myMap1 = new map<string, TH1F*>();
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
  
  smhtt_2017 t(argv[1],argv[2], isMC);
  t.Loop(maxEvents,reportEvery, SampleName , isMC);
  //delete myMap1;
  cout<<" Outpt written to "<<outputfile<<endl;
  sw.Stop();
  sw.Print();
  //delete myMap1;
  return 0;
}

void smhtt_2017::Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
{
   if (fChain == 0) return;
   double eventWeight=1.0;
   bool isMC = false;
   if( _isMC_=="MC" ) { isMC=true; }
   else  { isMC=false;  }
   
   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   //   ofstream myfile;
   // if(isMC)
   //   myfile.open ("eventAnalysis.txt");
   // myfile << "event"  << "\t" << "lumi" << "\t" << "run" <<
   //   "\t" << "ele pt" <<
   //   "\t" << "tau pt" <<
   //   "\t" << "met"    <<
   //   "\t" << "taudm"  <<
   //   "\t" << "tauGM"  <<
   //   "\n";
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      eventWeight=1.0;
      //cout<<"entry :"<<jentry<<endl;
      makeMyPlot("a", 0, 0, 0, eventWeight);
      bool trigger=((passDoubleTightTau35TightID && pt_1>40 && pt_2>40) 
		    || (passDoubleMediumTau40TightID && pt_1>40 && pt_2>40) 
		    || (passDoubleTightTau40ID && pt_1>40 && pt_2>40));
      //bool triggerHPS=(passDoubleMediumHPSTau35 && pt_1>40 && pt_2>40);
      //bool trigger2016=((passDoubleTau35 or passDoubleTauCmb35) && pt_1>40 && pt_2>40);

      //if(selectTrgEle27 || selectTrgEle32 || selectTrgEle35 || selectTrgCross)
      if( trigger )
	{
	  makeMyPlot("b", 0, 0, 0, eventWeight);
	  TLorentzVector mytau2; 
	  mytau2.SetPtEtaPhiM(pt_2,eta_2,phi_2,m_2);
	  TLorentzVector mytau1;
	  mytau1.SetPtEtaPhiM(pt_1,eta_1,phi_1,m_1);
	  
	  int mydm1=l1_decayMode;
	  int mydm2=l2_decayMode;
	  float gen1=gen_match_1;
	  float gen2=gen_match_2;
	  float signalRegion=(byMediumDeepVSjet_1);
	  float antiisoRegion=(byVVVLooseDeepVSjet_1 && !byMediumDeepVSjet_1);
	  float iso2=(byMediumDeepVSjet_2);
	  float iso1=(byMediumDeepVSjet_1);
	  	  

	  if(mytau1.Pt()>40 && fabs(mytau1.Eta())<2.1 && iso1>0.5
	     && byVVVLooseDeepVSe_1 && byVLooseDeepVSmu_1
	     )
	    {
	      makeMyPlot("c", 0, 0, 0, eventWeight);
	      // selct taus
	      if(mytau2.Pt()>40 && fabs(mytau2.Eta())<2.1 && iso2>0.5
		 && byVVVLooseDeepVSe_2 && byVLooseDeepVSmu_2 
		 )
		{
		  makeMyPlot("d", 0, 0, 0, eventWeight);
		  // opp charge
		  if(q_1*q_2 <0)
		    {
		      makeMyPlot("e", 0, 0, 0, eventWeight);
		      //cout<<"eventWeight="<<eventWeight<<endl;
		      double applySf=1.0;
		      // if(isMC)
		      // 	applySf=getScaleFactors( pt_1, pt_2, eta_1, eta_2, l2_decayMode);
		      
		      TLorentzVector dau1, dau2;
		      dau1=mytau1; dau2=mytau2;
		      double mt=TMass_F(dau1.Pt(), dau1.Phi() , met, metphi);
		      float mvis=(dau1+dau2).M();
		      eventWeight =  eventWeight*applySf;
		      // match filers
		      makeMyPlot("f", 0, 0, 0, eventWeight);
		      
		      
		      // 3rd lepton veto
		      makeMyPlot("g", 0, 0, 0, eventWeight);
		      
		      //bjet veto
		      //if(nbtag < 1)
		      //if(nbtag < 1 && nbtagL <2)
			{
			  makeMyPlot("h", 0, 0, 0, eventWeight);
			  
			  //dr cut
			  double dr_etau=dau1.DeltaR(dau2);
			  if(dr_etau>0.5)
			    {
			      makeMyPlot("i", 0, 0, 0, eventWeight);
			      
			      //mT cut
			      double mT_tauMet=TMass_F(dau2.Pt(), dau2.Phi() , met, metphi);
			      if(mT_tauMet<50)
				{
				  makeMyPlot("j", 0, 0, 0, eventWeight);
		
				}
			    }
			}
		      
		    }
		}
	    }
	}
   }
   //myfile.close();
  
}
void smhtt_2017::BookHistos(const char* file1, const char* file2)
{
  TFile* file_in= new TFile(file1, "READ");
  
  fileName = new TFile(file2, "RECREATE");
  fileName->cd();
  
}
double smhtt_2017::delta_R(float phi1, float eta1, float phi2, float eta2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar   = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}


float smhtt_2017::pTvecsum_F(TLorentzVector a, TLorentzVector b, TLorentzVector c) {
  float pt_vecSum = (a + b+ c).Pt();
  return pt_vecSum;
}
double smhtt_2017::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi  -= 2.0*pi;
  if(dphi<= -1*pi) dphi +=  2.0*pi;
  return fabs(dphi);
}

float smhtt_2017::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
  //return  sqrt(2.0*LepPt*met*(1.0-cos(DeltaPhi(LepPhi, metPhi))));
  return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

void smhtt_2017::makeMyPlot( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight){
  string hNumber = histNumber;
  TLorentzVector mytau2; 
  mytau2.SetPtEtaPhiM(pt_2,eta_2,phi_2,m_2);
  TLorentzVector mytau1;
  mytau1.SetPtEtaPhiM(pt_1,eta_1,phi_1,m_1);
  
  int mydm1=l1_decayMode;
  int mydm2=l2_decayMode;
  float gen1=gen_match_1;
  float gen2=gen_match_2;
  float signalRegion=(byMediumDeepVSjet_1);
  float antiisoRegion=(byVVVLooseDeepVSjet_1 && !byMediumDeepVSjet_1);
  float iso2=(byMediumDeepVSjet_2);
  
  
  plotFill("tau1Pt_"+hNumber,  mytau1.Pt() , 38 , 24 , 100,  event_weight);
  plotFill("tau1Eta_"+hNumber, mytau1.Eta(), 25, -2.5, 2.5,  event_weight);
  plotFill("tau1Phi_"+hNumber, mytau1.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("tau2Pt_"+hNumber,  mytau2.Pt() , 25 , 30 , 80 ,  event_weight);
  plotFill("tau2Eta_"+hNumber, mytau2.Eta(), 25, -2.5, 2.5,  event_weight);
  plotFill("tau2Phi_"+hNumber, mytau2.Phi(), 30, -3.14, 3.14,  event_weight);
  
  plotFill("tau1DecayMode_"+hNumber,  mydm1 ,  12, 0, 12, event_weight);
  plotFill("tau2DecayMode_"+hNumber,  mydm2 ,  12, 0, 12, event_weight);
  //decayModeFinding_2
  plotFill("decayModeFinding_"+hNumber, decayModeFinding_2 , 4, 0, 2,  event_weight);
  plotFill("nJet_"+hNumber, njets , 8, 0, 8,  event_weight);
  
  int triggerBin1, triggerBin2, triggerBin3, triggerBin4;
  triggerBin1=triggerBin2=triggerBin3=triggerBin4=0;
  if( passDoubleTightTau40ID && pt_1>40 && pt_2>40 )  triggerBin3=3;
  if( passDoubleMediumTau40TightID && pt_1>40 && pt_2>40 )  triggerBin2=2;
  if( passDoubleTightTau35TightID && pt_1>40 && pt_2>40 ) triggerBin1=1;
  if(triggerBin1>0)
    plotFill("trigger_"+hNumber, triggerBin1 , 5, 0, 5,  event_weight);
  if(triggerBin2>0)
    plotFill("trigger_"+hNumber, triggerBin2 , 5, 0, 5,  event_weight);
  if(triggerBin3>0)
    plotFill("trigger_"+hNumber, triggerBin3 , 5, 0, 5,  event_weight);
  if(triggerBin4>0)
    plotFill("trigger_"+hNumber, triggerBin4 , 5, 0, 5,  event_weight);
  
  plotFill("mass1_"+hNumber, m_1 , 50, 0, 0.05,  event_weight);
  plotFill("mass2_"+hNumber, m_2 , 25, 0, 2.5,  event_weight);
  plotFill("genmatch1_"+hNumber, gen_match_1 , 7, 0, 7,  event_weight);
  plotFill("genmatch2_"+hNumber, gen_match_2 , 7, 0, 7,  event_weight);
  
  double deltaPhi = DeltaPhi(mytau1.Phi() , mytau2.Phi());
  double deltaEta = fabs(mytau1.Eta() - mytau2.Eta());
  plotFill("deltaPhi_"+hNumber, deltaPhi , 30, -3.14, 3.14,  event_weight);
  plotFill("deltaEta_"+hNumber, deltaEta ,  25, -2.5, 2.5,  event_weight);
  double deltaR = mytau1.DeltaR(mytau2);
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);
  TLorentzVector corrMet;
  corrMet.SetPtEtaPhiE(met,0 ,metphi, met);
  double higgsPt = pTvecsum_F(mytau1, mytau2, corrMet);
  plotFill("higgsPt_"+hNumber, higgsPt , 25, 0, 250,  event_weight);
  plotFill("met_"+hNumber, met , 20, 0, 100,  event_weight);
  plotFill("metPhi_"+hNumber, metphi , 30, -3.14, 3.14,  event_weight);
  double mT_eleMet = TMass_F( mytau2.Pt(), mytau2.Phi(), met, metphi  );
  plotFill("mT_eMet_"+hNumber,  mT_eleMet ,  30, 0, 150,  event_weight);
  int nEvent=1;
  plotFill("nEvents_"+hNumber, nEvent , 3, 0.0, 3.0 , event_weight);
  plotFill("eventWeight_"+hNumber, event_weight , 20, 0.0 , 2.0,  1.0 );
  plotFill("event_"+hNumber, evt/1000000, 1000, 0.0, 1000,  event_weight);
  //cout<<"     elePt_"<<hNumber<<" = "<< elePt->at(tmpCand[0])<<endl;
  // TLorentzVector dau1, dau2;
  // dau1.SetPtEtaPhiE(pt_1, eta_1, phi_1, e_1);
  // dau2.SetPtEtaPhiE(pt_2, eta_2, phi_2, e_2);

}

double smhtt_2017::getScaleFactors( double elept, double taupt, double eleeta, double taueta, int taudm)
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
  double sf_taufesSF = 1.0;
  
  eleRecoSF_corr=h_eleRecoSF_highpt->GetBinContent(h_eleRecoSF_highpt->GetXaxis()->FindBin(eleeta),h_eleRecoSF_highpt->GetYaxis()->FindBin(elept));
  if (debug==true ) std::cout<<"eleRecoSF_corr =  "<< eleRecoSF_corr<<std::endl;
  eleEffSF_corr=h_eleIDSF->GetBinContent(h_eleIDSF->GetXaxis()->FindBin(eleeta),h_eleIDSF->GetYaxis()->FindBin(elept));
  if (debug==true ) std::cout<<"eleEffSF_corr =  "<< eleEffSF_corr<<std::endl;
  if (debug==true ) std::cout<<"This works line 269 "<<std::endl;
  

  if( gen_match_2>=5)
    {
      sf_tauidSF_m = h_tauidSF_m->GetBinContent(h_tauidSF_m->GetXaxis()->FindBin(taudm));
      //sf_tauidSF_m = fn_tauIDSF_m->Eval(tau_Pt->at(reco_tau));
      sf_tauidSF_vvvl = h_tauidSF_vvvl->GetBinContent(h_tauidSF_vvvl->GetXaxis()->FindBin(taudm));
      //sf_tauidSF_vvvl = fn_tauIDSF_vvl->Eval(tau_Pt->at(reco_tau));
      sf_tauesSF = h_tauesSF->GetBinContent(h_tauesSF->GetXaxis()->FindBin(taudm));
    }
  
  if(gen_match_2==1 || gen_match_2==3)
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
  if( gen_match_2==2 || gen_match_2==4 ){  ///  muons to pass deep tau 
    if(taudm==0)
      {
	if(abs(taueta) < 0.4 ) sf_fakeMu=1.14;
	if(abs(taueta) > 0.4 && abs(taueta) < 0.8 ) sf_fakeMu=1.0;
	if(abs(taueta) > 0.8 && abs(taueta) < 1.2 ) sf_fakeMu=0.87;
	if(abs(taueta) > 1.2 && abs(taueta) < 1.7 ) sf_fakeMu=0.52;
	if(abs(taueta) > 1.7 && abs(taueta) < 2.3 ) sf_fakeMu=1.47;
      }
    if(taudm==1)
      {
	if(abs(taueta) > 0.0 && abs(taueta) < 0.4 ) sf_fakeMu=0.69;
      }
  }

  TLorentzVector my_eleP4, my_tauP4, my_metP4;
  my_eleP4.SetPtEtaPhiE(pt_1, eta_1, phi_1 , e_1);
  my_tauP4.SetPtEtaPhiE(pt_2, eta_2, phi_2 , e_2);
  my_metP4.SetPtEtaPhiE(met , 0    , metphi, met);

  double higgsPt = pTvecsum_F(my_eleP4, my_tauP4, my_metP4);
  double higgPt_weight=1.0;
  if (njets==0)
    higgPt_weight = gr_NNLOPSratio_pt_mcatnlo_0jet->Eval(min(higgsPt,125.0));
  else if (njets==1)
    higgPt_weight = gr_NNLOPSratio_pt_mcatnlo_1jet->Eval(min(higgsPt,625.0));
  else if (njets==2)
    higgPt_weight = gr_NNLOPSratio_pt_mcatnlo_2jet->Eval(min(higgsPt,800.0));
  else if (njets>=3)
    higgPt_weight = gr_NNLOPSratio_pt_mcatnlo_3jet->Eval(min(higgsPt,925.0));
  else
    higgPt_weight = 1.0;
  
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
  // if(found_DYjet_sample)
  //   zptmass_weight= get_zptmass_weight();
  
  // double top_pt_weight=1.0;
  // if(found_TTbar_sample){
  //   int t_index = get_t_Cand(); int tbar_index = get_tbar_Cand();
  //   if( t_index >-1 && tbar_index > -1 ){
  //     top_pt_weight = sqrt( exponential(mcPt->at(t_index)) * exponential(mcPt->at(tbar_index)) );
  //     //cout<<"top_pt_weight = "<<top_pt_weight<<endl;
  //   } 
  // }
  
  
  if( elept<28.0 ) e_trg_sf=1.0;
  //sf_htt_workspace=sf_htt_workspace * e_trk_sf * e_idiso_sf * e_trg_sf * e_trg24_sf * t_trg_sf * t_deepid_tightvsele_sf;
  sf_htt_workspace=  e_trk_sf * e_idiso_sf *  e_trg24_sf * e_trg_sf;
  //  rv_sf = eleRecoSF_corr * eleEffSF_corr * sf_Zvtx * sf_tauidSF_m * sf_tauesSF *  sf_fakeMu * sf_taufesSF *sf_tauTrg* sf_htt_workspace;
  rv_sf = eleEffSF_corr * sf_tauidSF_m * sf_tauTrg * sf_fakeEle * sf_fakeMu  * sf_htt_workspace * higgPt_weight ;
  if(rv_sf>0)
    return rv_sf;
  else
    return 1.0;
  //return 1.0;
}
