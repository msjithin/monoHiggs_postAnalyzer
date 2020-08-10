
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
  t.Loop(maxEvents,reportEvery, SampleName , isMC);
  //delete myMap1;
  cout<<" Outpt written to "<<outputfile<<endl;
  sw.Stop();
  sw.Print();
  return 0;
}

void etau_analyzer::Loop(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
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
  
  int out_event_list[185] = {800, 4433346, 4413207, 4416347, 4519381, 12806227, 15617232, 16281684, 21209941, 21210321, 420580, 420588, 15510421, 1649067, 9175093, 867202, 13515877, 13517295, 823466, 824918, 15614219, 9126011, 9128537, 17259838, 9139876, 9173848, 9176415, 12669676, 12672835, 12674321, 16055569, 41801463, 41803065, 6006830, 70947887, 70835015, 70835251, 70837446, 429983, 440277, 12005634, 12005686, 12009397, 9662270, 464922, 466859, 468151, 472870, 481688, 490541, 28483758, 15934298, 16011815, 16014433, 42636353, 41597903, 41604340, 41604979, 15010857, 11030207, 11031503, 11032546, 7559039, 7570535, 14884510, 14897821, 729961, 11616947, 14919206, 455021, 814707, 7540353, 7540459, 640903, 2915496, 2834792, 2837510, 59694588, 53201, 5865475, 5868261, 22103609, 2890356, 617526, 58661287, 58661233, 59876261, 643779, 646082, 3007953, 58690637, 58700765, 58701244, 58705810, 503133, 1965929, 1967816, 1968240, 2017728, 2024942, 1975537, 1976450, 2028859, 9038280, 2035355, 2174509, 60093009, 60100542, 500192, 9917971, 25873325, 66793704, 66794119, 66794387, 66794709, 3869034, 22136033, 551314, 66778064, 66778514, 6630363, 6620027, 21755063, 66985928, 6675021, 6675303, 9013203, 9014134, 11315239, 13690187, 13693581, 6792536, 5832472, 19941066, 20047550, 20050080, 22156716, 22158950, 22091167, 24473453, 20051225, 20053470, 1554913, 5674609, 5674700, 20202950, 20334988, 612357, 613741, 32030654, 6500838, 6503948, 6508211, 32086851, 32088835, 6682240, 32097739, 6686167, 1038954, 5764535, 1207989, 1344229, 1404374, 19909345, 19412576, 22342666, 22343067, 22345224, 22101146, 767499, 13477345, 13479549, 13486971, 19410353, 29832676, 3193144, 3195098, 3199969, 19494818, 14945755, 14947562, 2820093, 112480341, 2947155, 2950608};

  int good_event_list[200] = { 12473, 4436039, 4438116, 4439801, 427301, 12806341, 395374, 17255299, 17255801, 16288238, 423352, 15516289, 15511055, 15537219, 15537474, 15548827, 1645188, 1647935, 12136693, 12654875, 15513313, 869765, 873828, 860701, 9136111, 9138304, 9173004, 20916, 12668470, 41797249, 41811942, 6004977, 6007124, 76412757, 70834828, 70841582, 4851877, 5998514, 438560, 475206, 12090067, 5316498, 15934842, 16014287, 42626022, 42636920, 15002074, 76413747, 38848, 43947, 7543260, 7543917, 7551385, 7552702, 11032715, 7559973, 7568135, 7570098, 11618517, 633266, 634640, 14865541, 14891947, 739968, 14899227, 14908724, 745235, 742479, 813909, 815432, 597822, 6019745, 2837149, 59692705, 59693099, 22103754, 2889937, 2892125, 58660739, 58663638, 59875959, 58697603, 501097, 504316, 1962431, 1979617, 2033933, 2034008, 2179110, 60093925, 76630069, 76631091, 66984015, 3866930, 66798895, 83591, 22137384, 22138278, 22140900, 22140985, 22148229, 8724310, 66775877, 66783957, 66786969, 66787537, 66788435, 66796027, 717675, 726778, 729472, 6618453, 9014022, 11314971, 13689994, 5834238, 19943934, 19944709, 19994819, 22340162, 22157306, 22094572, 16786036, 616219, 32029849, 32030877, 32031053, 6503365, 6506019, 6513657, 32087774, 32090909, 32091587, 32092046, 6647740, 1032251, 6684662, 1040136, 1042276, 1044482, 30126299, 4579208, 5766504, 16299077, 1407399, 22101308, 649516, 14597095, 773982, 788486, 3135955, 3141620, 19400402, 19409494, 19410842, 3191443, 3212753, 14939009, 112696104, 112702333, 2821688, 1628696, 22608879, 22613627, 69351665, 69407207, 113545065, 39360463, 112544017, 112544526, 112545580, 112552621, 112558014, 39449674, 112566848, 100682577, 100689551, 100699702, 79570991, 79848057, 79850998, 79851449, 1246273, 5909009, 37849262, 694969, 79662721, 79668722, 79670891, 1366898, 96608660, 5959521, 1435354, 7175094, 20310401, 699089, 703140, 20323316, 15921585};

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
   ofstream myfile;
   ofstream myfile2;
   ofstream myfile3;
   ofstream ptcd;
   //if(isMC)
   ptcd.open("compare_ptdiff_etau.txt");
   //ptcd << "ele pt" <<"\n";
   
   myfile.open ("eventAnalysisEtau_"+sample+".txt");
   myfile << "event" << std::setw(10) << "lumi" << std::setw(10) << "run" << 
     std::setw(10) << "ele pt" << std::setw(10) << "ele eta" << std::setw(10) << "ele phi" << 
     std::setw(10) << "tau pt" << std::setw(10) << "tau eta" << std::setw(10) << "tau phi" << 
     std::setw(10) << "nEle"   << std::setw(10) << "nTau"    << std::setw(10) << "nMu" <<
     std::setw(10) << "nMC"    << std::setw(10) << "genMatch2"<<   
     std::setw(10) << "diEleVeto" <<std::setw(10) << "eVetoZtt" << std::setw(10) << "mVetoZtt" << std::setw(10)<<     
     "\n";

   myfile2.open ("rogue_events_"+sample+".txt");
   myfile2 << "event" << std::setw(10) << "lumi" << std::setw(10) << "run" <<
     std::setw(10) << "ele pt" << std::setw(10) << "ele eta" << std::setw(10) << "ele phi" <<
     std::setw(10) << "tau pt" << std::setw(10) << "tau eta" << std::setw(10) << "tau phi" <<
     std::setw(10) << "nEle"   << std::setw(10) << "nTau"    << std::setw(10) << "nMu" << 
     std::setw(10) << "nMC"    << std::setw(10) << "genMatch2"<<   
     std::setw(10) << "diEleVeto" <<std::setw(10) << "eVetoZtt" << std::setw(10) << "mVetoZtt" << std::setw(10)<<
     "\n";
   
   myfile3.open ("good_events_"+sample+".txt");
   myfile3 << "event" << std::setw(10) << "lumi" << std::setw(10) << "run" <<
     std::setw(10) << "ele pt" << std::setw(10) << "ele eta" << std::setw(10) << "ele phi" <<
     std::setw(10) << "tau pt" << std::setw(10) << "tau eta" << std::setw(10) << "tau phi" <<
     std::setw(10) << "nEle"   << std::setw(10) << "nTau"    << std::setw(10) << "nMu" << 
     std::setw(10) << "nMC"    << std::setw(10) << "genMatch2"<<   
     std::setw(10) << "diEleVeto" <<std::setw(10) << "eVetoZtt" << std::setw(10) << "mVetoZtt" << std::setw(10)<<
     "\n";

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
       
       eleTau_selector=true;
       eleEle_selector=true;
       ////Trigger bit selection
       
       if( ( (HLTEleMuX>>3&1 == 1 )      //HLT_Ele27_WPTight_Gsf_v
	     || (HLTEleMuX>>61&1 == 1)  //HLT_Ele32_WPTight_Gsf_v
	     || (HLTEleMuX>>5&1 == 1)   //HLT_Ele35_WPTight_Gsf_v
	     ))
       	 passSingleTriggerPaths=true;  //
       
       if( ( HLTTau>>1&1 == 1 ) )      //HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1
       	 passCrossTrigger=true;
       

       
       /////
       if(debug)cout<<"entry # : "<<jentry<<endl;
       event_weight=weight;
       if(debug)cout<<"reco selections begin"<<endl;
       eleCand.clear(); ele2Cand.clear();  tauCand.clear();
       ////// reco selection begin
       
       bool Ztt_selector=false;
       if(isMC)
	 event_weight=2.447677534;
       else
	 event_weight=1.0;

       eleCand.clear(); ele2Cand.clear();  tauCand.clear();
       if(eleTau_selector)
	 {
	   
	   if(metFilters==0)
	     {
	       if(debug)cout<<"metfilters selected"<<endl;
	       if (isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	       nMETFiltersPassed+=event_weight;
	       makeMyPlot("a", -1,0,0,event_weight);
	       if(debug)cout<<"genweight applied"<<endl;
	       if( passSingleTriggerPaths || passCrossTrigger  )
		 {
		   nSingleTrgPassed+=event_weight;
		   if(debug)cout<<"trigger selected"<<endl;
		   makeMyPlot("b", -1,0,0,event_weight);
		   eleCand = getEleCand(24.0,2.1);  ///// ele selected 
		   if( eleCand.size() >0 ) 
		     { 
		       nGoodMuonPassed+=event_weight;
		       if(debug)cout<<"this worked Line 526"<<endl;
		       
		       makeMyPlot("c", eleCand[0],0,0,event_weight);
		       tauCand = getTauCand(30.0,2.3);
		       if( tauCand.size() >0 )
			 {
			   nGoodTauPassed+=event_weight;
			   
			   reco_ele.clear();reco_tau.clear();
			   reco_ele=eleCand; reco_tau=tauCand;
			   makeMyPlot("d", reco_ele[0],0,0,event_weight);
			   if(found_DYjet_sample){
			     if( passDiElectronVeto(eleCand[0])==true 
				 && eVetoZTTp001dxyz(eleCand[0], tauCand[0])
				 && mVetoZTTp001dxyz(eleCand[0], tauCand[0])
				 ) Ztt_selector=true;
			     else Ztt_selector=false;
			   }
			   else
			     Ztt_selector=true;
			   
			   if(Ztt_selector) 
			     {
			       if (  eleCharge->at(reco_ele[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
				 {
				   nGoodMuTauPassed+=event_weight;
				   
				   if(debug)cout<<"this worked Line 538"<<endl;
				   afterSF1+=event_weight;
				   makeMyPlot("e", reco_ele[0],0,0,event_weight);
				   if ( MatchTriggerFilter(reco_ele[0], reco_tau[0]) )
				     {
				       if(debug)cout<<"this worked Line 534"<<endl;
				       
				       // if(debug)cout<<" sf : "<<getScaleFactors( reco_ele[0] , reco_tau[0] , false , isMC , debug ) <<endl;
				       // if (isMC) event_weight = event_weight * getScaleFactors( reco_ele[0] , reco_tau[0] , false , isMC , debug );
				       
				       afterSF4+=event_weight;
				       makeMyPlot("f", reco_ele[0],0,0,event_weight);
				       if( thirdLeptonVeto(reco_ele[0] , reco_tau[0]) < 0 )
					 {
					   nPassedThirdLepVeto+=event_weight;
					   makeMyPlot("g", reco_ele[0],0,0,event_weight);
					   if( passBjetVeto(reco_ele[0] , reco_tau[0], -1) == true)
					     {
					       nPassedBjetVeto+=event_weight;
					       makeMyPlot("h", reco_ele[0],0,0,event_weight);
					       
					       double deltaR = dR(eleEta->at(reco_ele[0]), elePhi->at(reco_ele[0]),  tau_Eta->at(reco_tau[0]), tau_Phi->at(reco_tau[0]));
					       if(deltaR > 0.5 )
						 {
						   
						   
						   // for(int imc=0; imc<nMC; imc++){
						   //   if(abs(mcPID->at(imc))==15 || abs(mcPID->at(imc))==111 || abs(mcPID->at(imc))==211 )
						   //     cout<<"mcPID : "<<mcPID->at(imc)<<" , genMatch2 :"
						   // 	 << (genMatch2->at(imc)>>1&1==1) <<" "
						   // 	 << (genMatch2->at(imc)>>2&1==1)<<" "
						   // 	 << (genMatch2->at(imc)>>3&1==1)<<" "
						   // 	 << (genMatch2->at(imc)>>4&1==1)<<" "
						   // 	 << (genMatch2->at(imc)>>5&1==1)<<" "
						   // 	 << (genMatch2->at(imc)>>6&1==1)<<" "
						   // 	 <<endl; 
						   // }						
						   
						   bool found_tauEle=false; bool found_tauTau=false;
						   std::vector<int> eleGenCand;  eleGenCand.clear();
						   for(int iMC=0;iMC<nMC;iMC++) //Loop over mc
						     {
						       if(fabs(mcPID->at(iMC))==11 && mcStatus->at(iMC)==62) eleGenCand.push_back(iMC);
						     }
						   std::vector<int> tauGenCand;  tauGenCand.clear();
						   for(int iMC=0;iMC<nMC;iMC++) //Loop over mc
						     {
						       if(fabs(mcPID->at(iMC))==15 && mcStatus->at(iMC)==62) tauGenCand.push_back(iMC);
						     }
						   if(eleGenCand.size()>0 && tauGenCand.size()>0)
						     {
						       double mc_dr = delta_R(mcPhi->at(eleGenCand[0]), mcEta->at(eleGenCand[0]), mcPhi->at(tauGenCand[0]), mcEta->at(tauGenCand[0]) );
						       if(mc_dr<0.2)
							 found_tauEle=true;
						     }
						   if(tauGenCand.size()>1)
						     {
						       double mc_dr_tt = delta_R(mcPhi->at(tauGenCand[0]), mcEta->at(tauGenCand[0]), mcPhi->at(tauGenCand[1]), mcEta->at(tauGenCand[1]) );
						       if(mc_dr_tt<0.2)
							 found_tauTau=true;
						     }
						   //cout<<"event:"<<event<<"  genmatch3="<<found_tauEle<< "  genmatch5="<<found_tauTau<<endl;
						   std::vector<int> v = gen_matching();
						   int genMatch2_=-1;
						   if(v.size() >0)
						     {
						       genMatch2_=v[0];
						       // cout<<"event:"<<event<<"   v: ";
						       // for(int imc=0; imc<v.size(); imc++){
						       //   cout<<v[imc]<<" ";
						       // }
						       // cout<<""<<endl;
						     }
						   
						   int diEleVeto = 0; int eVetoZtt=0; int mVetoZtt=0;
						   if( passDiElectronVeto(reco_ele[0])==true ) diEleVeto=1;
						   if(  eVetoZTTp001dxyz(eleCand[0], tauCand[0]) ) eVetoZtt = 1;
						   if(   mVetoZTTp001dxyz(eleCand[0], tauCand[0])) mVetoZtt=1;
						   
						   
						   
						   myfile << event << std::setw(10) << lumis << std::setw(10) << run << 
						     std::setw(10) << elePt->at(reco_ele[0]) << std::setw(10) << eleEta->at(reco_ele[0]) << std::setw(10) << elePhi->at(reco_ele[0]) << 
						     std::setw(10) << tau_Pt->at(reco_tau[0]) << std::setw(10) << tau_Eta->at(reco_tau[0]) << std::setw(10) <<  tau_Phi->at(reco_tau[0]) <<
						     std::setw(10) << nEle   << std::setw(10) << nTau    << std::setw(10) << nMu <<
						     std::setw(10) << nMC    << std::setw(10) << genMatch2_ << std::setw(10) <<
						     diEleVeto  << std::setw(10) << eVetoZtt  << std::setw(10) << mVetoZtt  << std::setw(10) <<
						     "\n";
						   
						 if (std::find(std::begin(out_event_list), std::end(out_event_list), event) != std::end(out_event_list))
						   {
						     //cout<<"rougue event found "<<event<<endl;
						     myfile2 << event << std::setw(10) << lumis << std::setw(10) << run <<
						       std::setw(10) << elePt->at(reco_ele[0]) << std::setw(10) << eleEta->at(reco_ele[0]) << std::setw(10) << elePhi->at(reco_ele[0]) <<
						       std::setw(10) << tau_Pt->at(reco_tau[0]) << std::setw(10) << tau_Eta->at(reco_tau[0]) << std::setw(10) <<  tau_Phi->at(reco_tau[0]) <<
						       std::setw(10) << nEle   << std::setw(10) << nTau    << std::setw(10) << nMu <<
						       std::setw(10) << nMC    << std::setw(10) << genMatch2_ << std::setw(10) << 
						       diEleVeto  << std::setw(10) << eVetoZtt  << std::setw(10) << mVetoZtt  << std::setw(10) <<
						       "\n";
						   }
						 else
						   {
						     //cout<<"rougue event found "<<event<<endl;
						     myfile3 << event << std::setw(10) << lumis << std::setw(10) << run <<
						       std::setw(10) << elePt->at(reco_ele[0]) << std::setw(10) << eleEta->at(reco_ele[0]) << std::setw(10) << elePhi->at(reco_ele[0]) <<
						       std::setw(10) << tau_Pt->at(reco_tau[0]) << std::setw(10) << tau_Eta->at(reco_tau[0]) << std::setw(10) <<  tau_Phi->at(reco_tau[0]) <<
						       std::setw(10) << nEle   << std::setw(10) << nTau    << std::setw(10) << nMu <<
						       std::setw(10) << nMC    << std::setw(10) << genMatch2_ << std::setw(10) <<
						       diEleVeto  << std::setw(10) << eVetoZtt  << std::setw(10) << mVetoZtt  << std::setw(10) <<
						       "\n";
						   }
						 float relEleIso = ( elePFChIso->at(reco_ele[0]) + max( elePFNeuIso->at(reco_ele[0]) + elePFPhoIso->at(reco_ele[0]) - 0.5 *elePFPUIso->at(reco_ele[0]) , 0.0 )) / (elePt->at(reco_ele[0]));
						 int eleid=-1;
						 if(eleIDbit->at(reco_ele[0])>>8&1==1)
						   eleid=1;
						 if (std::find(std::begin(good_event_list), std::end(good_event_list), event) != std::end(good_event_list))
						   {
						     ptcd <<event<<"\t"<< elePt->at(reco_ele[0]) <<"\t"<< eleEta->at(reco_ele[0])<<"\t"<< elePhi->at(reco_ele[0])
							  << "\t"      << tau_Pt->at(reco_tau[0])<<"\t"<< tau_Eta->at(reco_tau[0])<<"\t"<< tau_Phi->at(reco_tau[0])<<"\t"
							  << eleid <<"\t" << (relEleIso) <<"\t"
							  << (tau_byMediumDeepTau2017v2p1VSjet->at(reco_tau[0])==1) <<"\t" << (tau_byTightDeepTau2017v2p1VSe->at(reco_tau[0])==1)  <<"\t" << (tau_byVLooseDeepTau2017v2p1VSmu->at(reco_tau[0])==1)<<"\t"
							  << eleCharge->at(reco_ele[0]) <<"\t" << tau_Charge->at(reco_tau[0]) <<"\t"
							  << deltaR  <<"\t"<< genMatch2_ <<  "\t"
							  <<"\n";
						     makeMyPlot("cd", reco_ele[0],0,0,event_weight);
						   }
						 nDeltaRPassed+=event_weight;
						 if(isMC==false)event_weight=1.0;
						 if(debug)cout<<"this worked Line 558"<<endl;
						 //fillHist("5", reco_ele[0], reco_tau[0], event_weight, isMC);
						 makeMyPlot("i", reco_ele[0],0,0,event_weight);
						 double mT_eleMet = TMass_F((elePt->at(reco_ele[0])),(elePhi->at(reco_ele[0])),pfMET,pfMETPhi  );
						 if( mT_eleMet < 50 )
						   {
						     //fillHist("6", reco_ele[0], reco_tau[0], event_weight, isMC);
						     makeMyPlot("j", reco_ele[0],0,0,event_weight);
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
   cout<<"dr="<<delta_R(-1.2972, -0.765188,  0.515507, 0.350012)<< "  expected: 1.1152"<<endl;
   myfile.close();
   myfile2.close();
   myfile3.close();
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
   
   std::cout<<"HLTTau>>1&1 ==1 = "<<hltele61<<std::endl;
   std::cout<<"eleF45 = "<<eleF45<<std::endl;
   std::cout<<"eleF53 = "<<eleF53<<std::endl;
   std::cout<<"eleF54 = "<<eleF54<<std::endl;
   std::cout<<"eleF55 = "<<eleF55<<std::endl;
   std::cout<<"tauF11 = "<<tauF11<<std::endl;
   std::cout<<"tauF12 = "<<tauF12<<std::endl;
   std::cout<<"tauF16 = "<<tauF16<<std::endl;

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
  //Loop over electrons                                                                     
  for(int iEle=0;iEle<nEle;iEle++)
    {
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
	  //&& tau_LeadChargedHadron_dz->at(iTau) < 0.2
	  && fabs(tau_Charge->at(iTau))==1
  	  )kinematic = true;
      if(  tau_byVVVLooseDeepTau2017v2p1VSjet->at(iTau)==1 && tau_byMediumDeepTau2017v2p1VSjet->at(iTau)!=1 ) tauIsolation=true;
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
std::vector<int> etau_analyzer::getZCand()
{
  std::vector<int> tmpCand;  tmpCand.clear();
  for(int iMC=0;iMC<nMC;iMC++) //Loop over mc
    {
      if(fabs(mcPID->at(iMC))==23 && mcStatus->at(iMC)==62) tmpCand.push_back(iMC); 
    }
  return tmpCand;
}
std::vector<int> etau_analyzer::getJetCand(int eleIndex, int tauIndex, int ele2Index){
  std::vector<int> tmpCand;  tmpCand.clear();
  for(int iJet=0;iJet<nJet;iJet++) //Loop over jets
    {
      bool kinematic30 = false;
      bool kinematic50 = false; bool kinematic50Loose = false;
      bool jet_id = false; bool drPassed=false;
      if( jetPt->at(iJet) > 50 
	  && abs(jetEta->at(iJet))<2.65
	  && abs(jetEta->at(iJet))>3.139
	  && (jetID->at(iJet)>>0&1)==1
	  ) kinematic50=true;
      // else if( jetPt->at(iJet) < 50  
      // 	       && jetPUFullID->at(iJet)>>1&1==1 
      // 	       ) kinematic50Loose=true;
      else if( jetPt->at(iJet) > 30
	       && abs(jetEta->at(iJet))<4.7
	       && (jetID->at(iJet)>>0&1)==1
	       ) kinematic30=true;
      double lepton1Phi=elePhi->at(eleIndex);
      double lepton1Eta= eleEta->at(eleIndex);
      double lepton2Phi=0;double lepton2Eta=0;
      if(ele2Index<0){ lepton2Phi= tau_Phi->at(tauIndex); lepton2Eta=tau_Eta->at(tauIndex); }
      else           { lepton2Phi= elePhi->at(ele2Index); lepton2Eta=eleEta->at(ele2Index); }
      
      double dr_jetEle=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton1Phi, lepton1Eta );
      double dr_jetTau=delta_R( jetPhi->at(iJet), jetEta->at(iJet) , lepton2Phi, lepton2Eta);
      if( dr_jetEle>0.5 && dr_jetTau>0.5 )
	drPassed=true;
	  
      if( (kinematic50 || kinematic30 ) && drPassed==true && jetPUFullID->at(iJet)>>1&1==1)
	tmpCand.push_back(iJet);
    }
  return tmpCand;
}
//The noisy jets are defined as: 20 < pt < 50 && abs(eta) > 2.65 && abs(eta) < 3.139. 
bool etau_analyzer::noisyJet2017(){

  bool noisyJet = false;
  for(int iJet=0;iJet<nJet;iJet++) //Loop over jets
    {
      bool kinematic = false;
      if( jetPt->at(iJet) > 20 
	  && jetPt->at(iJet) < 50
	  && abs(jetEta->at(iJet)) > 2.65
          && abs(jetEta->at(iJet)) < 3.139
	  ) kinematic=true;

      if( kinematic )
        noisyJet=true;
    }
  return noisyJet;
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
      //if( muIDbit->at(iMu)>>6&1==1) relative_iso =true;
      if(muonId==true && kinematic==true && relative_iso==true){
	tmpCand.push_back(iMu);
      }                   
    }          
  if(tmpCand.size() > 0){ thirdLepIndex = tmpCand[0]; thirdLepVeto=false;}
  return thirdLepIndex;
  
}
int etau_analyzer::thirdLeptonVeto(int eleIndex, int tauIndex){
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
	found_3rdmu=true;
      }
    }

  if(tmpCand.size() > 0 && found_3rdmu){ thirdLepIndex = tmpCand[0]; thirdLepVeto=false;}
  return thirdLepIndex;
  
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
float etau_analyzer::pTvecsum_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float pt_vecSum = (a + b+ met).Pt();
  return pt_vecSum;
}

bool etau_analyzer::passBjetVeto(int eleIndex, int tauIndex, int ele2Index)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  double lepton1Phi=elePhi->at(eleIndex);
  double lepton1Eta= eleEta->at(eleIndex);
  double lepton2Phi=0;double lepton2Eta=0;
  if(ele2Index<0){ lepton2Phi= tau_Phi->at(tauIndex); lepton2Eta=tau_Eta->at(tauIndex); }
  else           { lepton2Phi= elePhi->at(ele2Index); lepton2Eta=eleEta->at(ele2Index); }
  
  double dr_jetEle=0.0; double dr_jetTau=0.0;
  
  for(int iJets=0; iJets<nJet ; iJets++){
    dr_jetEle=delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton1Phi, lepton1Eta );
    dr_jetTau=delta_R( jetPhi->at(iJets), jetEta->at(iJets) , lepton2Phi, lepton2Eta);
    
    if( jetPt->at(iJets) > 25  
	&& abs(jetEta->at(iJets)) < 2.4 
	&& jetID->at(iJets)>>0&1==1 
	&& jetPUFullID->at(iJets)>>1&1==1
	&& dr_jetEle>0.5 && dr_jetTau>0.5
	){
      tmpJetCand.push_back(iJets);
    }
  }
  //cout<<" jet size = "<<tmpJetCand.size()<<endl;
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
    if (deltaR > 0.15 && eleCharge->at(iElePlus[0])*eleCharge->at(iEleMinus[0])<0) 
      {
	return false;
      }
  }
  return true;
  
}
bool etau_analyzer::eVetoZTTp001dxyz(int eleIndex, int tauIndex){
  std::vector<int> tmpCand;
  tmpCand.clear();
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
      double deltaR_et = delta_R(tau_Phi->at(tauIndex), tau_Eta->at(tauIndex), elePhi->at(tmpCand[0]), eleEta->at(tmpCand[0]));
      double deltaR_ee = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), elePhi->at(tmpCand[0]), eleEta->at(tmpCand[0]));
      if(deltaR_et>0.5 && deltaR_ee>0.5)
	return false;
    }
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

std::vector<int> etau_analyzer::gen_matching(){
  int tmpCand=-1;
  std::vector<int> tmpGenMatch;
  tmpGenMatch.clear();
  
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>1&1==1 ) tmpGenMatch.push_back(1);
    else if( genMatch2->at(imc)>>2&1==1 ) tmpGenMatch.push_back(2);
    else if( genMatch2->at(imc)>>3&1==1 ) tmpGenMatch.push_back(3);
    else if( genMatch2->at(imc)>>4&1==1 ) tmpGenMatch.push_back(4);
    else if( genMatch2->at(imc)>>5&1==1 ) tmpGenMatch.push_back(5);
    else if( genMatch2->at(imc)>>6&1==1 ) tmpGenMatch.push_back(6);
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
double etau_analyzer::getScaleFactors( int eleIndex , int tauIndex, bool fakeBkg , bool isMC, bool debug)
{
  return 1.0;
}
bool etau_analyzer::MatchTriggerFilter(int eleIndex, int tauIndex)
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool passFilter = true;
  bool eleTriggerFilterMatch=false;
  int nEleTriggerFilterMatch=0;
  bool tauTriggerFilterMatch=false;
  int nTauTriggerFilterMatch=0;
  int nEleDoubleTriggerFilterMatch=0; bool eleDoubleTrgsMatch=false;
  for(int ifilter=39;ifilter<56;ifilter++)
    {
      if(eleFiredSingleTrgs->at(eleIndex)>>ifilter&1==1)
	{
	  eleTriggerFilterMatch=true;
	  nEleTriggerFilterMatch++;
	}
    }
  for(int ifilter=12;ifilter<20;ifilter++)
    {
      if(eleFiredDoubleTrgs->at(eleIndex)>>ifilter&1==1)
        {
          eleDoubleTrgsMatch=true;
          nEleDoubleTriggerFilterMatch++;
        }
    }
  if( 
     (HLTEleMuX>>3&1 == 1 && eleFiredSingleTrgs->at(eleIndex)>>12&1==1)
     || (HLTEleMuX>>61&1 == 1  && (eleFiredSingleTrgs->at(eleIndex)>>51&1==1 || eleFiredSingleTrgs->at(eleIndex)>>52&1==1) )
     || (HLTEleMuX>>5&1 == 1  && eleFiredSingleTrgs->at(eleIndex)>>14&1==1)
     || (HLTTau>>1&1 ==1  )//&& (eleFiredSingleTrgs->at(eleIndex)>>53&1==1 || eleFiredSingleTrgs->at(eleIndex)>>54&1==1 || eleFiredSingleTrgs->at(eleIndex)>>55&1==1))
      )
    return true;
  else
    return false;
}

double  etau_analyzer::getFR(int tauIndex){
    double frWeight=1.0;
  double tau_FR = 1.0;
  double tauPt=0.0;
  if( tau_Pt->at(tauIndex) < 120 )
    tauPt=tau_Pt->at(tauIndex);
  else
    tauPt=119.0;
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

void etau_analyzer::fillHist( string histNumber , int eleIndex, int tauIndex, float event_weight, bool isMC){
  string hNumber = histNumber;
  plotFill("elePt_"+hNumber,  elePt->at(eleIndex) , 38 , 24 , 100,  event_weight);
  plotFill("eleEta_"+hNumber, eleEta->at(eleIndex), 48, -2.4, 2.4,  event_weight);
  plotFill("elePhi_"+hNumber, elePhi->at(eleIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("eleDz_"+hNumber,  eleDz->at(eleIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("eleD0_"+hNumber,  eleD0->at(eleIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("electronID_"+hNumber, eleIDbit->at(eleIndex)>>8&1, 4, -2, 2,  event_weight); // electronID
  float relEleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (elePt->at(eleIndex));
  plotFill("relEleIso_"+hNumber, relEleIso, 15, 0, 0.3,  event_weight);
  plotFill("eleCharge_"+hNumber, eleCharge->at(eleIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 45, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), tau_Phi->at(tauIndex), tau_Eta->at(tauIndex));
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(eleIndex, tauIndex, -1);
  plotFill("nJet_"+hNumber,  jetCand.size() , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, pfMET , 20, 0, 100,  event_weight);
  
  double mT_eleMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
  plotFill("mT_eleMet_"+hNumber, mT_eleMet , 30, 0, 150,  event_weight);

  TLorentzVector myTau; 
  myTau.SetPtEtaPhiE(tau_Pt->at(tauIndex),tau_Eta->at(tauIndex),tau_Phi->at(tauIndex), tau_Energy->at(tauIndex));
  TLorentzVector myEle; 
  myEle.SetPtEtaPhiE(elePt->at(eleIndex), eleEta->at(eleIndex), elePhi->at(eleIndex), eleE->at(eleIndex));
  double visMass_mutau = VisMass_F(myTau, myEle);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  double HiggsPt = pTvecsum_F(myEle, myTau, myMet );
  plotFill("higgsPt_"+hNumber,HiggsPt , 25, 0, 250,  event_weight);

  double tot_tr_mass = TotTMass_F(myEle, myTau, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  int triggerBin=0;
  if( HLTTau>>1&1 == 1 )     triggerBin=4;    //HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1
  else if( HLTEleMuX>>3&1 == 1 ) triggerBin=1;   //HLT_Ele27_WPTight_Gsf_v
  else if( HLTEleMuX>>61&1 == 1 ) triggerBin=2;   //HLT_Ele32_WPTight_Gsf_v
  else if( HLTEleMuX>>5&1 == 1 ) triggerBin=3;     //HLT_Ele35_WPTight_Gsf_v
  plotFill("trigger_"+hNumber, triggerBin , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(isMC){
    if(found_GenMatch(1)) genMatchBin=1;
    else if(found_GenMatch(2)) genMatchBin=2;
    else if(found_GenMatch(3)) genMatchBin=3;
    else if(found_GenMatch(4)) genMatchBin=4;
    else if(found_GenMatch(5)) genMatchBin=5;
    else if(found_GenMatch(6)) genMatchBin=6;
  }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}
void etau_analyzer::fillHist( string histNumber , TLorentzVector eleP4, TLorentzVector tauP4, int eleIndex, int tauIndex, float event_weight, bool isMC){
  string hNumber = histNumber;
  plotFill("elePt_"+hNumber,  elePt->at(eleIndex) , 38 , 24 , 100,  event_weight);
  plotFill("eleEta_"+hNumber, eleEta->at(eleIndex), 48, -2.4, 2.4,  event_weight);
  plotFill("elePhi_"+hNumber, elePhi->at(eleIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("eleDz_"+hNumber,  eleDz->at(eleIndex), 40, -0.5, 0.5,  event_weight);
  plotFill("eleD0_"+hNumber,  eleD0->at(eleIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("electronID_"+hNumber, eleIDbit->at(eleIndex)>>8&1, 4, -2, 2,  event_weight); 
  float relEleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (elePt->at(eleIndex));
  plotFill("relEleIso_"+hNumber, relEleIso, 15, 0, 0.3,  event_weight);
  plotFill("eleCharge_"+hNumber, eleCharge->at(eleIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 45, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), tau_Phi->at(tauIndex), tau_Eta->at(tauIndex));
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);

  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(eleIndex, tauIndex, -1);
  plotFill("nJet_"+hNumber, jetCand.size() , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, pfMET , 20, 0, 100,  event_weight);

  double mT_eleMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
  plotFill("mT_eleMet_"+hNumber, mT_eleMet , 30, 0, 150,  event_weight);
  
  double visMass_mutau = VisMass_F(eleP4, tauP4);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  double HiggsPt = pTvecsum_F(eleP4, tauP4, myMet);
  plotFill("higgsPt_"+hNumber,HiggsPt , 25, 0, 250,  event_weight);

  double tot_tr_mass = TotTMass_F(eleP4, tauP4, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  int triggerBin=0;
  if( HLTTau>>1&1 == 1 )     triggerBin=4;
  else if( HLTEleMuX>>3&1 == 1 ) triggerBin=1;
  else if( HLTEleMuX>>61&1 == 1 ) triggerBin=2;
  else if( HLTEleMuX>>5&1 == 1 ) triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(isMC){
    if(found_GenMatch(1)) genMatchBin=1;
    else if(found_GenMatch(2)) genMatchBin=2;
    else if(found_GenMatch(3)) genMatchBin=3;
    else if(found_GenMatch(4)) genMatchBin=4;
    else if(found_GenMatch(5)) genMatchBin=5;
    else if(found_GenMatch(6)) genMatchBin=6;
  }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
}
void etau_analyzer::fillHist_dyll( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight, bool isMC){
  string hNumber = histNumber;
  plotFill("elePt_"+hNumber,  elePt->at(eleIndex) , 38 , 24 , 100,  event_weight);
  plotFill("eleEta_"+hNumber, eleEta->at(eleIndex), 48, -2.4, 2.4,  event_weight);
  plotFill("elePhi_"+hNumber, elePhi->at(eleIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("eleDz_"+hNumber,  eleDz->at(eleIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("eleD0_"+hNumber,  eleD0->at(eleIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("electronID_"+hNumber, eleIDbit->at(eleIndex)>>8&1, 4, -2, 2,  event_weight); // electronID
  float relEleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (elePt->at(eleIndex));
  plotFill("relEleIso_"+hNumber, relEleIso, 15, 0, 0.3,  event_weight);
  plotFill("eleCharge_"+hNumber, eleCharge->at(eleIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 45, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  
  double deltaR = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), elePhi->at(ele2Index), eleEta->at(ele2Index));
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(eleIndex, -1, ele2Index);
  plotFill("nJet_"+hNumber,  jetCand.size() , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, pfMET , 20, 0, 100,  event_weight);
  
  double mT_eleMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
  plotFill("mT_eleMet_"+hNumber, mT_eleMet , 30, 0, 150,  event_weight);

  TLorentzVector myTau; 
  myTau.SetPtEtaPhiE(tau_Pt->at(tauIndex),tau_Eta->at(tauIndex),tau_Phi->at(tauIndex), tau_Energy->at(tauIndex));
  TLorentzVector myEle; 
  myEle.SetPtEtaPhiE(elePt->at(eleIndex), eleEta->at(eleIndex), elePhi->at(eleIndex), eleE->at(eleIndex));
  TLorentzVector myEle2;
  myEle2.SetPtEtaPhiE(elePt->at(ele2Index), eleEta->at(ele2Index), elePhi->at(ele2Index), eleE->at(ele2Index));
  double visMass_eletau = VisMass_F(myEle2, myEle);
  plotFill("visMass_"+hNumber, visMass_eletau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  double HiggsPt = pTvecsum_F(myEle, myEle2, myMet );
  plotFill("higgsPt_"+hNumber,HiggsPt , 25, 0, 250,  event_weight);

  double tot_tr_mass = TotTMass_F(myEle, myEle2, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  int triggerBin=0;
  if( HLTTau>>1&1 == 1 )     triggerBin=4;
  else if( HLTEleMuX>>3&1 == 1 ) triggerBin=1;
  else if( HLTEleMuX>>61&1 == 1 ) triggerBin=2;
  else if( HLTEleMuX>>5&1 == 1 ) triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(isMC)
    {
      if(found_GenMatch(1)) genMatchBin=1;
      else if(found_GenMatch(2)) genMatchBin=2;
      else if(found_GenMatch(3)) genMatchBin=3;
      else if(found_GenMatch(4)) genMatchBin=4;
      else if(found_GenMatch(5)) genMatchBin=5;
      else if(found_GenMatch(6)) genMatchBin=6;
    }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);  
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
void etau_analyzer::makeMyPlot( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight){
  string hNumber = histNumber;
  int ele_select=-1;
  std::vector<int> tmpCand; tmpCand.clear();
  for(int iEle=0;iEle<nEle;iEle++)
    {
      tmpCand.push_back(iEle);
    }
  if(eleIndex>-1)
    ele_select=eleIndex;
  else
    ele_select=tmpCand[0];
  plotFill("elePt_"+hNumber,  elePt->at(ele_select) , 38 , 24 , 100,  event_weight);
  //cout<<"     elePt_"<<hNumber<<" = "<< elePt->at(tmpCand[0])<<endl;
}
