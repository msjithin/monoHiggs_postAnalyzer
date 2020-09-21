#define smhet_2017_cxx
#include "smhet_2017.h"
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
  
  smhet_2017 t(argv[1],argv[2], isMC);
  t.Loop(maxEvents,reportEvery, SampleName , isMC);
  //delete myMap1;
  cout<<" Outpt written to "<<outputfile<<endl;
  sw.Stop();
  sw.Print();
  //delete myMap1;
  return 0;
}

void smhet_2017::Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
{
   if (fChain == 0) return;
   double eventWeight=1.0;
   bool isMC = false;
   if( _isMC_=="MC" ) { isMC=true; }
   else  { isMC=false;  }
   
   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   ofstream myfile;
   ofstream ptcd;
   int good_event_list[200] = { 12473, 4436039, 4438116, 4439801, 427301, 12806341, 395374, 17255299, 17255801, 16288238, 423352, 15516289, 15511055, 15537219, 15537474, 15548827, 1645188, 1647935, 12136693, 12654875, 15513313, 869765, 873828, 860701, 9136111, 9138304, 9173004, 20916, 12668470, 41797249, 41811942, 6004977, 6007124, 76412757, 70834828, 70841582, 4851877, 5998514, 438560, 475206, 12090067, 5316498, 15934842, 16014287, 42626022, 42636920, 15002074, 76413747, 38848, 43947, 7543260, 7543917, 7551385, 7552702, 11032715, 7559973, 7568135, 7570098, 11618517, 633266, 634640, 14865541, 14891947, 739968, 14899227, 14908724, 745235, 742479, 813909, 815432, 597822, 6019745, 2837149, 59692705, 59693099, 22103754, 2889937, 2892125, 58660739, 58663638, 59875959, 58697603, 501097, 504316, 1962431, 1979617, 2033933, 2034008, 2179110, 60093925, 76630069, 76631091, 66984015, 3866930, 66798895, 83591, 22137384, 22138278, 22140900, 22140985, 22148229, 8724310, 66775877, 66783957, 66786969, 66787537, 66788435, 66796027, 717675, 726778, 729472, 6618453, 9014022, 11314971, 13689994, 5834238, 19943934, 19944709, 19994819, 22340162, 22157306, 22094572, 16786036, 616219, 32029849, 32030877, 32031053, 6503365, 6506019, 6513657, 32087774, 32090909, 32091587, 32092046, 6647740, 1032251, 6684662, 1040136, 1042276, 1044482, 30126299, 4579208, 5766504, 16299077, 1407399, 22101308, 649516, 14597095, 773982, 788486, 3135955, 3141620, 19400402, 19409494, 19410842, 3191443, 3212753, 14939009, 112696104, 112702333, 2821688, 1628696, 22608879, 22613627, 69351665, 69407207, 113545065, 39360463, 112544017, 112544526, 112545580, 112552621, 112558014, 39449674, 112566848, 100682577, 100689551, 100699702, 79570991, 79848057, 79850998, 79851449, 1246273, 5909009, 37849262, 694969, 79662721, 79668722, 79670891, 1366898, 96608660, 5959521, 1435354, 7175094, 20310401, 699089, 703140, 20323316, 15921585};

   
   // if(isMC)
   //   myfile.open ("eventAnalysis.txt");
   // myfile << "event" << "\t" << "lumi" << "\t" << "run" << 
   //   "\t" << "ele pt" << "\t" << "ele eta" << "\t" << "ele phi" << 
   //   "\t" << "tau pt" << "\t" << "tau eta" << "\t" << "tau phi" << 
   //    "\t" << "nEle"   << "\t" << "nTau"    << "\t" << "nMu" <<
   //   "\t" << "nMC"    << "\t" << "genMatch2"<<   "\t" <<
   //   "diEleVeto"    << "\t" << "eVetoZtt"<<   "\t" << "mVetoZtt" <<   "\t" <<
   //   "\n";
   // if(isMC)
   //   ptcd.open("compare_ptdiff_smhtt.txt");
   // //ptcd << "ele pt" <<"\n";
   // //eventWeight=2.447677534; // to normalize for cross section
   // // if(isMC)
   // //   eventWeight=2.447677534;
   // // else
   // //eventWeight=1.0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      eventWeight=1.0;
      //cout<<"entry :"<<jentry<<endl;
      makeMyPlot("a", 0, 0, 0, eventWeight);
      bool selectTrgEle27=(passEle27 && matchEle27_1 && filterEle27_1 && pt_1>28);
      bool selectTrgEle32=(passEle32 && matchEle32_1 && filterEle32_1 && pt_1>33);
      bool selectTrgEle35=(passEle35 && matchEle35_1 && filterEle35_1 && pt_1>36);
      bool selectTrgCross=(passEle24Tau30 && matchEle24Tau30_1 && filterEle24Tau30_1
      			   && matchEle24Tau30_2 && filterEle24Tau30_2
      			   && pt_1>25 && pt_1<28 && pt_2>35 && abs(eta_1)<2.1 && abs(eta_2)<2.1
      			   );
      // bool selectTrgEle27=(passEle27 && pt_1>28 && abs(eta_1)<2.1);
      // bool selectTrgEle32=(passEle32 && pt_1>33 && abs(eta_1)<2.1);
      // bool selectTrgEle35=(passEle35 && pt_1>36 && abs(eta_1)<2.1);
      // bool selectTrgCross=(passEle24Tau30 && pt_1>25 && pt_1<28 && pt_2>35  
      // 			   && abs(eta_1)<2.1 && abs(eta_2)<2.1);
      if(selectTrgEle27 || selectTrgEle32 || selectTrgEle35 || selectTrgCross)
      //if( selectTrgEle32 || selectTrgEle35 || selectTrgCross)
	{
	  makeMyPlot("b", 0, 0, 0, eventWeight);
	  
	  //select electrons
	  if(pt_1>25 && fabs(eta_1)<2.1 && iso_1<0.15)// && eid90_noiso_1 )
	    {
	      makeMyPlot("c", 0, 0, 0, eventWeight);
	      // selct taus
	      if(pt_2>30 && fabs(eta_1)<2.3 && byMediumDeepVSjet_2>0
		 && byTightDeepVSe_2>0.5 && byVLooseDeepVSmu_2>0.5 
		 && (l2_decayMode==0 || l2_decayMode==1 || l2_decayMode==10 || l2_decayMode==11)
		 //&& decayModeFinding_2>0.5
		 )
		{
		  makeMyPlot("d", 0, 0, 0, eventWeight);
		  // opp charge
		  if(q_1*q_2 <0)
		    {
		      makeMyPlot("e", 0, 0, 0, eventWeight);
		      //cout<<"eventWeight="<<eventWeight<<endl;
		      double applySf=1.0;
		      if(isMC)
			applySf=getScaleFactors( pt_1, pt_2, eta_1, eta_2, l2_decayMode);
		      
		      eventWeight =  eventWeight*applySf;
		      // match filers
		      makeMyPlot("f", 0, 0, 0, eventWeight);
		      
		      
		      // 3rd lepton veto
		      makeMyPlot("g", 0, 0, 0, eventWeight);
		      
		      //bjet veto
		      if(nbtag <= 0 && nbtagL <= 1)
			{
			  makeMyPlot("h", 0, 0, 0, eventWeight);
			  
			  //dr cut
			  double dr_etau=dR(eta_1, phi_1, eta_2, phi_2);
			  if(dr_etau>0.5)
			    {
			      makeMyPlot("i", 0, 0, 0, eventWeight);
			      myfile << evt << "\t" << lumi << "\t" << run << 
				"\t" << pt_1 << "\t" << eta_1 << "\t" << phi_1 << 
                                "\t" << pt_2 << "\t" << eta_2 << "\t" << phi_2 <<
				"\t" << "nEle" << "\t" << "nTau" << "\t" <<"nMu" <<
				"\t" << dr_etau<<
				"\t" << "nMC" << "\t" << gen_match_2 << "\t"<<gen_match_1 << "\t"<<
				"\n";
			      int eleid=-1;
			      if(eid90_noiso_1)
				eleid=1;
			      if (std::find(std::begin(good_event_list), std::end(good_event_list), evt) != std::end(good_event_list))
				{
				  // ptcd << evt <<"\t" << pt_1 <<"\t" << eta_1 <<"\t" << phi_1 << "\t"
				  //      << pt_2 <<"\t" << eta_2 <<"\t" << phi_2  << "\t"
				  //      << eleid <<"\t" << iso_1<<"\t" 
				  //      << byMediumDeepVSjet_2 <<"\t" << byTightDeepVSe_2  <<"\t" << byVLooseDeepVSmu_2<<"\t"
				  //      << q_1 <<"\t" << q_2 <<"\t"
				  //      <<dr_etau  <<"\t"<< gen_match_2 <<  "\t" << gen_match_1 <<  "\t"
				  //      <<passEle27 <<"\t"<< passEle32 <<"\t"<< passEle35 <<"\t"<< passEle24Tau30 <<"\t"
				  //      <<matchEle27_1 <<"\t"<<matchEle32_1 <<"\t"<<matchEle35_1 <<"\t"<<matchEle24Tau30_1 <<"\t"<<matchEle24Tau30_2 <<"\t"
				  //      << filterEle27_1<<"\t"<< filterEle32_1<<"\t"<< filterEle35_1 <<"\t"<<filterEle24Tau30_1<<"\t"<< filterEle24Tau30_2<<"\t"
				  //      <<"\n";
				}
			      //mT cut
			      double mT_eMet=TMass_F(pt_1, phi_1 , met, metphi);
			      if(mT_eMet<50)
				{
				  makeMyPlot("j", 0, 0, 0, eventWeight);
				  //cout<<"run:"<<run<<" lumi:"<<lumi<<" event:"<<evt<<" genpX:"<<genpX<<" genpY:"<<genpY<<" vispX:"<<vispX<<" vispY"<<vispY<<endl;
				  //cout<<"l2_decayMode="<<l2_decayMode<<"  , decayModeFinding_2="<<decayModeFinding_2<<endl;
				}
			    }
			}
		      
		    }
		}
	    }
	}
   }
   myfile.close();
}
void smhet_2017::BookHistos(const char* file1, const char* file2)
{
  TFile* file_in= new TFile(file1, "READ");
  
  fileName = new TFile(file2, "RECREATE");
  fileName->cd();
  
}
double smhet_2017::delta_R(float phi1, float eta1, float phi2, float eta2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar   = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}



double smhet_2017::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
}

float smhet_2017::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
  //return  sqrt(2.0*LepPt*met*(1.0-cos(DeltaPhi(LepPhi, metPhi))));
  return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

void smhet_2017::makeMyPlot( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight){
  string hNumber = histNumber;
  plotFill("elePt_"+hNumber,  pt_1 , 38 , 24 , 100,  event_weight);
  plotFill("tauDecayMode_"+hNumber,  l2_decayMode ,  12, 0, 12, event_weight);
  //decayModeFinding_2
  plotFill("decayModeFinding_"+hNumber, decayModeFinding_2 , 4, 0, 2,  event_weight);
  plotFill("nJet_"+hNumber, njets , 6, 0, 6,  event_weight);
  
  int triggerBin=0;
  if( passEle35 && matchEle35_1 && filterEle35_1 ) triggerBin=4;
  else if( passEle32 && matchEle32_1 && filterEle32_1 ) triggerBin=3;
  else if( passEle27 && matchEle27_1 && filterEle27_1 ) triggerBin=2;
  else if( passEle24Tau30 && matchEle24Tau30_1 && matchEle24Tau30_2 && filterEle24Tau30_1 && filterEle24Tau30_2 ) triggerBin=1;
  plotFill("mass1_"+hNumber, m_1 , 20, 0, 0.2,  event_weight);
  plotFill("mass2_"+hNumber, m_2 , 25, 0, 2.5,  event_weight);
  plotFill("genmatch1_"+hNumber, gen_match_1 , 7, 0, 7,  event_weight);
  plotFill("genmatch2_"+hNumber, gen_match_2 , 7, 0, 7,  event_weight);

  plotFill("trigger_"+hNumber, triggerBin , 5, 0, 5,  event_weight);
  plotFill("met_"+hNumber, met , 20, 0, 100,  event_weight);
  plotFill("metPhi_"+hNumber, metphi , 30, -3.14, 3.14,  event_weight);
  double mT_eleMet = TMass_F( pt_1, phi_1, met, metphi  );
  plotFill("mT_eMet_"+hNumber,  mT_eleMet ,  30, 0, 150,  event_weight);

  //cout<<"     elePt_"<<hNumber<<" = "<< elePt->at(tmpCand[0])<<endl;
}

double smhet_2017::getScaleFactors( double elept, double taupt, double eleeta, double taueta, int taudm)
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
  

  if( gen_match_2==5)
    {
      sf_tauidSF_m = h_tauidSF_m->GetBinContent(h_tauidSF_m->GetXaxis()->FindBin(taudm));
      //sf_tauidSF_m = fn_tauIDSF_m->Eval(tau_Pt->at(reco_tau));
      sf_tauidSF_vvvl = h_tauidSF_vvvl->GetBinContent(h_tauidSF_vvvl->GetXaxis()->FindBin(taudm));
      //sf_tauidSF_vvvl = fn_tauIDSF_vvl->Eval(tau_Pt->at(reco_tau));
      sf_tauesSF = h_tauesSF->GetBinContent(h_tauesSF->GetXaxis()->FindBin(taudm));
    }
  
  if(gen_match_2<5)
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
  if(taudm==0 && abs(taueta)<=1.4 ) sf_taufesSF = h_taufesSF->Eval(1);
  if(taudm==0 && abs(taueta)>1.4 )  sf_taufesSF = h_taufesSF->Eval(3);
  if(taudm==1 && abs(taueta)<=1.4 ) sf_taufesSF = h_taufesSF->Eval(5);
  if(taudm==1 && abs(taueta)>1.4 )  sf_taufesSF = h_taufesSF->Eval(7);
  
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
  // if(genZCand.size()>0){
  //   w->var("z_gen_pt")->setVal(mcPt->at(genZCand[0]));
  //   w->var("z_gen_mass")->setVal(mcMass->at(genZCand[0]));
  //   //cout<< "zptmass_weight_nom = "<<w->function("zptmass_weight_nom")->getVal();
  //   zptmass_weight=w->function("zptmass_weight_nom")->getVal();
  // }

  if(elept<28)
    e_trg_sf=1.0;
  if(elept>28)
    e_trg24_sf=1.0;
  // if(taupt>35)
  //   t_trg_sf=1.0;
  //sf_htt_workspace=sf_htt_workspace * e_trk_sf * e_idiso_sf * e_trg_sf * e_trg24_sf * t_trg_sf * t_deepid_tightvsele_sf;
  sf_htt_workspace=  e_trk_sf * e_idiso_sf *  e_trg24_sf * e_trg_sf;
  //cout<<"sf_htt_workspace : "<<sf_htt_workspace<<endl;
  
  rv_sf = eleRecoSF_corr * eleEffSF_corr * sf_Zvtx * sf_tauidSF_m * sf_tauesSF *  sf_fakeMu * sf_taufesSF *sf_tauTrg* sf_htt_workspace;
  //rv_sf = sf_htt_workspace;
  //if(rv_sf < 0.5)
  //cout<<"rv_sf< 0.5 "<<" "<<eleRecoSF_corr<<" "<<eleEffSF_corr<<" "<<sf_Zvtx<<" "<<sf_tauidSF_vvvl<<" "<<sf_tauesSF<<" "<<sf_fakeMu<<" "<<sf_taufesSF<<" "<<sf_htt_workspace<<endl;
  if(rv_sf>0)
    return rv_sf;
  else
    return 1.0;
  //return 1.0;
}
