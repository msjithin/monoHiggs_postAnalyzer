#include <TH2.h>
//#include "ComputeWG1Unc.h"

#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include "TMultiGraph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <stdio.h>
#include <TF1.h>
#include <TDirectoryFile.h>
#include <TRandom3.h>
#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TKey.h"
#include "THashList.h"
#include "THStack.h"
#include "TPaveLabel.h"
#include "TFile.h"
//#include "myHelper.h"
//#include "tr_Tree.h"
//#include "ScaleFactor.h"
//#include "ZmmSF.h"
//#include "LumiReweightingStandAlone.h"
//#include "btagSF.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"
#include "TGraph2D.h"
#include "TColor.h"
#include <vector>
typedef std::vector<double> NumV;

using namespace std;
Double_t Luminosity = 41520.0;//Lumi for inclusive
bool debugOn=false;

/////        add the name of histograms to be created     
std::vector<string>  plotList = {
  "cutflow_n","cutflow_n_fr", "cutflow_n_dyll",
  "muPt_5", "muEta_5", "muPhi_5", "muDz_5", "muD0_5", "muonID_5", "relMuIso_5", "muCharge_5", "tauPt_5", "tauEta_5", "tauPhi_5", "tauIso_5", "tauDecayMode_5", "tauCharge_5", "tauAntiEle_5", "tauAntiMu_5" , "deltaR_5", "higgsPt_5", "nJet_5",  "visMass_5", "mT_muMet_5",  "met_5",  "tot_TMass_5", "trigger_5",
  
  "muPt_5_fr", "muEta_5_fr", "muPhi_5_fr", "muDz_5_fr", "muD0_5_fr", "muonID_5_fr", "relMuIso_5_fr", "muCharge_5_fr", "tauPt_5_fr", "tauEta_5_fr", "tauPhi_5_fr", "tauIso_5_fr", "tauDecayMode_5_fr", "tauCharge_5_fr", "tauAntiEle_5_fr", "tauAntiMu_5_fr" ,  "deltaR_5_fr", "higgsPt_5_fr", "nJet_5_fr", "visMass_5_fr", "mT_muMet_5_fr", "met_5_fr", "tot_TMass_5_fr", "trigger_5_fr",

  "muPt_6", "muEta_6", "muPhi_6", "muDz_6", "muD0_6", "muonID_6", "relMuIso_6", "muCharge_6", "tauPt_6", "tauEta_6", "tauPhi_6", "tauIso_6", "tauDecayMode_6", "tauCharge_6", "tauAntiEle_6", "tauAntiMu_6" , "deltaR_6", "higgsPt_6", "nJet_6",  "visMass_6", "mT_muMet_6",  "met_6",  "tot_TMass_6", "trigger_6",

  "muPt_6_fr", "muEta_6_fr", "muPhi_6_fr", "muDz_6_fr", "muD0_6_fr", "muonID_6_fr", "relMuIso_6_fr", "muCharge_6_fr", "tauPt_6_fr", "tauEta_6_fr", "tauPhi_6_fr", "tauIso_6_fr", "tauDecayMode_6_fr", "tauCharge_6_fr", "tauAntiEle_6_fr", "tauAntiMu_6_fr" ,  "deltaR_6_fr", "higgsPt_6_fr", "nJet_6_fr", "visMass_6_fr", "mT_muMet_6_fr", "met_6_fr", "tot_TMass_6_fr", "trigger_6_fr",

  "muPt_5_dyll", "muEta_5_dyll", "muPhi_5_dyll", "muDz_5_dyll", "muD0_5_dyll", "muonID_5_dyll", "relMuIso_5_dyll", "muCharge_5_dyll", "tauPt_5_dyll", "tauEta_5_dyll", "tauPhi_5_dyll", "tauIso_5_dyll", "tauDecayMode_5_dyll", "tauCharge_5_dyll", "tauAntiEle_5_dyll", "tauAntiMu_5_dyll" ,  "deltaR_5_dyll", "higgsPt_5_dyll", "nJet_5_dyll", "visMass_5_dyll", "mT_muMet_5_dyll", "met_5_dyll", "tot_TMass_5_dyll", "trigger_5_dyll",
  "muPt_6_dyll", "muEta_6_dyll", "muPhi_6_dyll", "muDz_6_dyll", "muD0_6_dyll", "muonID_6_dyll", "relMuIso_6_dyll", "muCharge_6_dyll", "tauPt_6_dyll", "tauEta_6_dyll", "tauPhi_6_dyll", "tauIso_6_dyll", "tauDecayMode_6_dyll", "tauCharge_6_dyll", "tauAntiEle_6_dyll", "tauAntiMu_6_dyll" ,  "deltaR_6_dyll", "higgsPt_6_dyll", "nJet_6_dyll", "visMass_6_dyll", "mT_muMet_6_dyll", "met_6_dyll", "tot_TMass_6_dyll", "trigger_6_dyll",

  "muPt_5_dyll_fr", "muEta_5_dyll_fr", "muPhi_5_dyll_fr", "muDz_5_dyll_fr", "muD0_5_dyll_fr", "muonID_5_dyll_fr", "relMuIso_5_dyll_fr", "muCharge_5_dyll_fr", "tauPt_5_dyll_fr", "tauEta_5_dyll_fr", "tauPhi_5_dyll_fr", "tauIso_5_dyll_fr", "tauDecayMode_5_dyll_fr", "tauCharge_5_dyll_fr", "tauAntiEle_5_dyll_fr", "tauAntiMu_5_dyll_fr" ,  "deltaR_5_dyll_fr", "higgsPt_5_dyll_fr", "nJet_5_dyll_fr", "visMass_5_dyll_fr", "mT_muMet_5_dyll_fr",  "met_5_dyll_fr", "tot_TMass_5_dyll_fr", "trigger_5_dyll_fr",
  "muPt_6_dyll_fr", "muEta_6_dyll_fr", "muPhi_6_dyll_fr", "muDz_6_dyll_fr", "muD0_6_dyll_fr", "muonID_6_dyll_fr", "relMuIso_6_dyll_fr", "muCharge_6_dyll_fr", "tauPt_6_dyll_fr", "tauEta_6_dyll_fr", "tauPhi_6_dyll_fr", "tauIso_6_dyll_fr", "tauDecayMode_6_dyll_fr", "tauCharge_6_dyll_fr", "tauAntiEle_6_dyll_fr", "tauAntiMu_6_dyll_fr" ,  "deltaR_6_dyll_fr", "higgsPt_6_dyll_fr", "nJet_6_dyll_fr", "visMass_6_dyll_fr", "mT_muMet_6_dyll_fr",   "met_6_dyll_fr", "tot_TMass_6_dyll_fr", "trigger_6_dyll_fr"

}; 

void make_hist(string input_file, string output_file, string histSaveName, string histname_string, string sample_name, TString directory_name ,  Double_t weight_lumi , bool isNLO)
{
  if(debugOn==true)cout<<"starting make_hist"<<endl;
  TString histname = TString(histname_string);
  TString sample = TString(sample_name); 
  TString input_name = TString(input_file);
  TString output_name = TString(output_file); 
  TString histSave_name = TString(histSaveName);

  Double_t scale_factor = 1.0;
  //***************** Tau id scale factor *****************
  
  Double_t net_weight = weight_lumi;
  
  if(sample=="data_obs"){ net_weight = 1.0;  }
  if(debugOn==true)cout<<"sample = "<< sample<<"  histname ="<<histname<<endl;
  TFile *file_input = new TFile(TString(input_name));
  TFile* outputFile;
  if(debugOn==true)cout<<"1. updating file for "<< histname<<endl;

  outputFile = new TFile(output_name, "UPDATE");
  TDirectory *dir_update = outputFile->mkdir(directory_name); 
  outputFile->cd(directory_name);
  bool histExists = file_input->GetListOfKeys()->Contains(histname);
  if(histExists==true){ 
    TH1F* histo_mutau = (TH1F*)((TH1F*)file_input->Get(histname))->Clone(TString(histSaveName+"_"+histname));
    if(sample_name=="data_obs"){ net_weight = 1.0;  }
    histo_mutau->Scale(net_weight);
    histo_mutau->Write();
    
  }

  file_input->Close();
  outputFile->Close();
  if(debugOn==true)cout<<".....................updated "<< histname<<endl;
}


int main(int argc, char** argv)
{

  if(debugOn==true)cout<<"starting main fn......."<<'\n';

  std::string input = *(argv + 1);
  std::string output = *(argv + 2);
  std::string sample = *(argv + 3);
  std::string histSaveName = *(argv + 4);
  
  float tes=0;
  if (argc > 1) {
    tes = atof(argv[5]);
  }
  TFile *f_Double = new TFile(input.c_str());
  cout<<left<<" Input file "<<std::setw(30)<<input.c_str()<<std::setw(10)<<" ---> "<<" Output written to  "<<std::setw(10)<<output.c_str()<<endl;
  //cout<<"  "<<endl;
  
  TH1F* nbevt;
  bool nEventsExists = f_Double->GetListOfKeys()->Contains("nEvents");
  if(nEventsExists){
    nbevt = (TH1F*) f_Double->Get("nEvents" );
  }
  else{
    cout<<"nEvents not found"<<endl;
    return 0;
  }

  //TTree *arbre = (TTree*) f_Double->Get("/ggNtuplizer/EventTree");
  //float ngen = 10000;
  if(debugOn==true)cout<<"Opened nEvents histogram "<<'\n';

  Double_t ngen = nbevt->GetBinContent(1);
  //  Double_t nFinal = nbevt->GetBinContent(12);

  if(debugOn==true)cout<<" "<<"ngen = "<< ngen<<"  "<<'\n';
  if (ngen<1) if(debugOn==true)cout<<"XXXXXXXXXXXXX  "<<" check this one !!" << input.c_str()<<"  XXXXXXXXXXXXX "<<'\n';

  Double_t xs=1.0; Double_t weight=1.0; Double_t luminosity=Luminosity;
  Double_t LOtoNNLO_DY = 5765.4/4954.0;
  Double_t LOtoNNLO_Wjets = 61526.7/50380; 
  if (sample=="ZL" || sample=="ZTT" || sample=="ZJ" || sample=="ZLL"){ xs=5765.4; weight=luminosity*xs/ngen;}
  else if (sample=="DY_LO" || sample=="DYJetsToLL"){ xs=LOtoNNLO_DY*4954.0; weight=luminosity*xs/ngen;}
  else if (sample=="DY1JetsToLL"){ xs=LOtoNNLO_DY*1012.5; weight=luminosity*xs/ngen;}
  else if (sample=="DY2JetsToLL"){ xs=LOtoNNLO_DY*332.8; weight=luminosity*xs/ngen;}
  else if (sample=="DY3JetsToLL"){ xs=LOtoNNLO_DY*101.8; weight=luminosity*xs/ngen;}
  else if (sample=="DY4JetsToLL"){ xs=LOtoNNLO_DY*54.8; weight=luminosity*xs/ngen;}
  else if (sample=="TTJets") {xs=831.76; weight=luminosity*xs/ngen;}
  else if (sample=="TTTo2L2Nu" ) {xs=88.29; weight=luminosity*xs/ngen;}
  else if (sample=="TTToHadronic" ) {xs=377.96; weight=luminosity*xs/ngen;}
  else if (sample=="TTToSemiLeptonic" ) {xs=365.35; weight=luminosity*xs/ngen;}

  else if (sample=="ZTTjet_inc"){ xs=LOtoNNLO_DY*4954.0; weight=2.447677534;}
  else if (sample=="ZTT1jet"){ xs=LOtoNNLO_DY*1012.5; weight=0.785055866;}
  else if (sample=="ZTT2jet"){ xs=LOtoNNLO_DY*332.8; weight=0.968954563;}
  else if (sample=="ZTT3jet"){ xs=LOtoNNLO_DY*101.8; weight=0.552080344;}
  else if (sample=="ZTT4jet"){ xs=LOtoNNLO_DY*54,8; weight=0.48914453;}
  
  else if (sample=="WJetsToLNu_inc") {xs=LOtoNNLO_Wjets*50380.0;    weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu")     {xs=LOtoNNLO_Wjets*50380.0;    weight=luminosity*xs/ngen;}
  else if (sample=="W1JetsToLNu")    {xs=LOtoNNLO_Wjets*9644.5;     weight=luminosity*xs/ngen;}
  else if (sample=="W2JetsToLNu")    {xs=LOtoNNLO_Wjets*3144.5;     weight=luminosity*xs/ngen;}
  else if (sample=="W3JetsToLNu")    {xs=LOtoNNLO_Wjets*964.8;      weight=luminosity*xs/ngen;}
  else if (sample=="W4JetsToLNu")    {xs=LOtoNNLO_Wjets*485.6;      weight=luminosity*xs/ngen;}

  else if (sample=="WJets_inc"){xs=LOtoNNLO_Wjets*50380.0;    weight=57.29217;}
  else if (sample=="W1Jet")    {xs=LOtoNNLO_Wjets*9644.5;     weight=7.715;}
  else if (sample=="W2Jet")    {xs=LOtoNNLO_Wjets*3144.5;     weight=17.0511;}
  else if (sample=="W3Jet")    {xs=LOtoNNLO_Wjets*964.8;      weight=2.357;}
  else if (sample=="W4Jet")    {xs=LOtoNNLO_Wjets*485.6;      weight=2.1368;}

  else if (sample=="WJetsToLNu_2J" ){xs=LOtoNNLO_Wjets*50380; weight=luminosity*xs/ngen;} // {xs=50380.0; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT100To200") {xs=LOtoNNLO_Wjets*1345.0; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT200To400") {xs=LOtoNNLO_Wjets*359.7; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT400To600") {xs=LOtoNNLO_Wjets*48.91; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT600To800") {xs=LOtoNNLO_Wjets*12.04; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT800To1200") {xs=LOtoNNLO_Wjets*5.52; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT1200To2500") {xs=LOtoNNLO_Wjets*1.33; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT2500ToInf") {xs=LOtoNNLO_Wjets*0.0322; weight=luminosity*xs/ngen;}

  else if (sample=="GluGluHToTauTau") {xs=48.58*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="VBFHToTauTau") {xs=3.782*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="ZHToTauTau") {xs=0.7612*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="WplusH125" || sample=="WplusHToTauTau") {xs=0.8400*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="WminusH125" || sample=="WminusHToTauTau") {xs=0.5328*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="GluGluHToWWTo2L2Nu") {xs=1.154; weight=luminosity*xs/ngen;}
  else if (sample=="VBFHToWWTo2L2Nu") {xs=0.0897; weight=luminosity*xs/ngen;}

  else if (sample=="GluGluZH_HToWW") {xs=0.0262; weight=luminosity*xs/ngen;}
  else if (sample=="HWminusJ_HToWW") {xs=0.114; weight=luminosity*xs/ngen;}
  else if (sample=="HWplusJ_HToWW") {xs=0.180; weight=luminosity*xs/ngen;}
  else if (sample=="HZJ_HToWW") {xs=0.163; weight=luminosity*xs/ngen;}
  else if (sample=="WGToLNuG") {xs=464.4; weight=luminosity*xs/ngen;}
  else if (sample=="ggZH_HToTauTau_ZToLL") {xs=0.1227*0.0627*3*0.033658; weight=luminosity*xs/ngen;}
  else if (sample=="ggZH_HToTauTau_ZToNuNu") {xs=0.1227*0.0627*0.2000; weight=luminosity*xs/ngen;}
  else if (sample=="ggZH_HToTauTau_ZToQQ") {xs=0.1227*0.0627*0.6991; weight=luminosity*xs/ngen;}
  else if (sample=="ttHToNonbb") {xs=1.0; weight=luminosity*xs/ngen;}

  else if (sample=="QCD") {xs=720648000*0.00042; weight=luminosity*xs/ngen;}
  else if (sample=="data_obs"){weight=1.0;}
 
  else if (sample=="ZJetsToNuNu_HT100-200") {xs=280.92; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT200-400") {xs=77.64; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT400-600") {xs=10.671; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT600-800") {xs=2.5611; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT800-1200") {xs=1.1778; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT1200-2500") {xs=0.2874; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT2500-Inf") {xs=0.006933; weight=luminosity*xs/ngen;}
  
  else if (sample=="WZTo1L3Nu") {xs=3.05; weight=luminosity*xs/ngen;}
  else if (sample=="WZTo1L1Nu2Q") {xs=10.71; weight=luminosity*xs/ngen;}
  else if (sample=="WZTo2L2Q") {xs=5.595; weight=luminosity*xs/ngen;}
  else if (sample=="WZTo3LNu") {xs=4.43; weight=luminosity*xs/ngen;}
  else if (sample=="VVTo2L2Nu") {xs=13.84; weight=luminosity*xs/ngen;}

  else if (sample=="ST_tW_antitop") {xs=35.85; weight=luminosity*xs/ngen;}
  else if (sample=="ST_tW_top") {xs=35.85; weight=luminosity*xs/ngen;}
  else if (sample=="ST_t-channel_antitop") {xs=80.95; weight=luminosity*xs/ngen;}
  else if (sample=="ST_t-channel_top") {xs=136.02; weight=luminosity*xs/ngen;}
  
  else if (sample=="WWTo1L1Nu2Q" || sample=="WWToLNuQQ") {xs=10.71; weight=luminosity*xs/ngen;}
  else if (sample=="WWTo2L2Nu") {xs=12.178; weight=luminosity*xs/ngen;}
  else if (sample=="GluGluWWTo2L2Nu") {xs=0.59; weight=luminosity*xs/ngen;}
  else if (sample=="WWTo2L2Nu_DoubleScattering") {xs=1.62; weight=luminosity*xs/ngen;} //WWTo2L2Nu_DoubleScattering
  else if (sample=="WpWpJJ_EWK_QCD") {xs=0.02615; weight=luminosity*xs/ngen;}
  else if (sample=="WpWpJJ_EWK") {xs=0.02615; weight=luminosity*xs/ngen;}
  else if (sample=="WpWpJJ_QCD") {xs=0.02615; weight=luminosity*xs/ngen;}

  else if (sample=="ZZTo2L2Q") {xs=3.22; weight=luminosity*xs/ngen;}
  else if (sample=="ZZTo2Q2Nu") {xs=4.03; weight=luminosity*xs/ngen;}
  else if (sample=="ZZTo4L") {xs=1.212; weight=luminosity*xs/ngen;}
  else if (sample=="ZZTo2L2Nu") {xs=0.564; weight=luminosity*xs/ngen;}
  
  else if (sample=="WWW") {xs=0.2086; weight=luminosity*xs/ngen;}
  else if (sample=="WWZ") {xs=0.1651; weight=luminosity*xs/ngen;}
  else if (sample=="WZZ") {xs=0.05565; weight=luminosity*xs/ngen;}
  else if (sample=="ZZZ") {xs=0.01398; weight=luminosity*xs/ngen;}

  else if (sample=="EWKWMinus" || sample=="EWKWMinus2Jets") {xs=23.24; weight=luminosity*xs/ngen;}
  else if (sample=="EWKWPlus" || sample=="EWKWPlus2Jets") {xs=29.59; weight=luminosity*xs/ngen;}
  else if (sample=="EWKZLL" || sample=="EWKZ2Jets_ZToLL" || sample=="EWKZLL_TT" || sample=="EWKZLL_J" || sample=="EWKZLL_L" || sample=="EWKZLL_LL") {xs=4.321; weight=luminosity*xs/ngen;}
  else if (sample=="EWKZNuNu" || sample=="EWKZ2Jets_ZToNuNu" || sample=="EWKZNuNu_TT" || sample=="EWKZNuNu_J" || sample=="EWKZNuNu_L" || sample=="EWKZNuNu_LL") {xs=10.66; weight=luminosity*xs/ngen;}

  else {
    cout<<"Attention!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    cout<<"***********                                  *****************"<<endl;
    cout<<"***********"<<sample.c_str()<<"  NOT found in sample list ************"<<endl;
    cout<<"***********                                  *****************"<<endl;
    cout<<"Attention!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    return 0;
  }
  
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(10);
  
  
  if(debugOn==true)cout<<"************************ works till here, after label *********************"<<endl;
  TString output_ = TString(output);
  TFile* f_output;
  f_output = new TFile(output_ , "RECREATE");
  std::vector<string> histnames;       histnames.clear();
  std::vector<TString> dirNames;      dirNames.clear();
  for(int i = 0; i < plotList.size(); i++){
    histnames.push_back(plotList[i]); //
    dirNames.push_back(plotList[i]);
    //TDirectory *dir_update = f_output->mkdir(plotList[i]);
  }
  if(debugOn==true)cout<<"This works too P2" << endl;
  f_output->Close();
  f_Double->Close();
  // if(debugOn==true)cout<<"This workks too before makehist" <<endl;
  //make_hist(input, output, "muPt_5", sample, dirName, weight, true, ngen , xs);
  for(int i = 0; i < histnames.size(); i++){
    make_hist(input, output, histSaveName, histnames[i] , sample, dirNames[i], weight, true );
  }
  if(debugOn==true)cout<<"Finished writing to "<<output << endl;
  
}



