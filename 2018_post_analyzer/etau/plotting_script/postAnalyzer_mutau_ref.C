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
Double_t Luminosity = 41521.0;//Lumi for inclusive
//Double_t Luminosity = 4823.0;    // for run B
//Double_t Luminosity = 9664.0;  // for run C
//Double_t Luminosity = 4252.0;  // for run D
//Double_t Luminosity = 9278.0;  // for run E
//Double_t Luminosity = 13540.0; // for run F

void make_hist(string input_file, string output_file, string histname_string, string sample_name, string name_ ,  Double_t weight_lumi , bool isNLO, double ngen_, double xsec)
{

  TString histname = TString(histname_string);
  TString sample = TString(sample_name); 
  TString input_name = TString(input_file);
  TString output_name = TString(output_file); 
  TString name = TString(name_); 
  // Double_t Luminosity = 44980.0;

  Double_t scale_factor = 1.0;
  //***************** Tau id scale factor *****************

  Double_t net_weight = weight_lumi;
      /*  if (sample=="WJetsToLNu" || sample=="WJetsToLNu_HT100To200" ||  sample=="WJetsToLNu_HT200To400" || sample=="WJetsToLNu_HT400To600" || sample=="WJetsToLNu_HT600To800" || sample=="WJetsToLNu_HT800To1200" || sample=="WJetsToLNu_HT1200To2500" || sample=="WJetsToLNu_HT2500ToInf" )
    {
      net_weight=net_weight*(1/0.89);
    }
  */
  if(sample=="data_obs"){ net_weight = 1.0;  }

  TFile *file_input = new TFile(TString(input_name));
  TFile* outputFile;
  // cout<<"1. updating file for "<< histname<<endl;

  outputFile = new TFile(output_name, "UPDATE");
  outputFile->cd(histname);
  
  if(histname == "Events_level_"){
    TH1F* histo_events_level = (TH1F*)((TH1F*)file_input->Get("Events_level"))->Clone(TString(name+"_"+histname));
    // cout<<"3. updating file for "<< histname<<endl;
    /*    if(sample=="data_obs"){
      histo_events_level->SetBinContent(10,0.0);
      histo_events_level->SetBinContent(11,0.0);
      histo_events_level->SetBinContent(12,0.0);
      }*/
    histo_events_level->Scale(net_weight);
    histo_events_level->Write();
  }
  else if(histname == "Cutflow_"){
    TH1F* histo_events_level = (TH1F*)((TH1F*)file_input->Get("Cutflow"))->Clone(TString(name+"_"+histname));
    // cout<<"3. updating file for "<< histname<<endl;
    histo_events_level->Scale(net_weight);
    histo_events_level->Write();
  }

  else if(histname == "MET_0_"){
    TH1F* histo_met_0 = (TH1F*)((TH1F*)file_input->Get("MET_0"))->Clone(TString(name+"_"+"MET_0"));
    TH1F* histo_met_1 = (TH1F*)((TH1F*)file_input->Get("MET_1"))->Clone(TString(name+"_"+"MET_1"));
    TH1F* histo_met_2 = (TH1F*)((TH1F*)file_input->Get("MET_2"))->Clone(TString(name+"_"+"MET_2"));
    TH1F* histo_met_3 = (TH1F*)((TH1F*)file_input->Get("MET_3"))->Clone(TString(name+"_"+"MET_3"));
    TH1F* histo_met_4 = (TH1F*)((TH1F*)file_input->Get("MET_4"))->Clone(TString(name+"_"+"MET_4"));

    // cout<<" updating file for "<< "MET "<<endl;
    histo_met_0->Scale(net_weight);
    histo_met_1->Scale(net_weight);
    histo_met_2->Scale(net_weight);
    histo_met_3->Scale(net_weight);
    histo_met_4->Scale(net_weight);

    histo_met_0->Write();
    histo_met_1->Write();
    histo_met_2->Write();
    histo_met_3->Write();
    histo_met_4->Write();

  }
  else{
    TH1F* histo_etau_0 = (TH1F*)((TH1F*)file_input->Get(histname+"0"))->Clone(TString(name+"_"+histname+"0"));
    TH1F* histo_etau_1 = (TH1F*)((TH1F*)file_input->Get(histname+"1"))->Clone(TString(name+"_"+histname+"1"));
    TH1F* histo_etau_2 = (TH1F*)((TH1F*)file_input->Get(histname+"2"))->Clone(TString(name+"_"+histname+"2"));
    TH1F* histo_etau_3 = (TH1F*)((TH1F*)file_input->Get(histname+"3"))->Clone(TString(name+"_"+histname+"3"));
    TH1F* histo_etau_4 = (TH1F*)((TH1F*)file_input->Get(histname+"4"))->Clone(TString(name+"_"+histname+"4"));
    TH1F* histo_etau_5 = (TH1F*)((TH1F*)file_input->Get(histname+"5"))->Clone(TString(name+"_"+histname+"5"));
    TH1F* histo_etau_6 = (TH1F*)((TH1F*)file_input->Get(histname+"6"))->Clone(TString(name+"_"+histname+"6"));
    //cout<<"declared histograms"<<endl;
    // cout<<"2. updating file for "<< histname<<endl;
    if(sample=="data_obs"){ net_weight = 1.0;  }
    if(sample=="data_obs" && histname =="pfMET_"){ net_weight = 0.0;  }
    histo_etau_0->Scale(net_weight);
    histo_etau_1->Scale(net_weight);
    histo_etau_2->Scale(net_weight);
    histo_etau_3->Scale(net_weight);
    histo_etau_4->Scale(net_weight);
    histo_etau_5->Scale(net_weight);
    histo_etau_6->Scale(net_weight);

    histo_etau_0->Write();
    histo_etau_1->Write();
    histo_etau_2->Write();
    histo_etau_3->Write();
    histo_etau_4->Write();
    histo_etau_5->Write();
    histo_etau_6->Write();
    
  }
  /* const int Nbins_0 = histo_etau_0->GetXaxis()->GetNbins();
    const int Nbins_1 = histo_etau_1->GetXaxis()->GetNbins();
    const int Nbins_2 = histo_etau_2->GetXaxis()->GetNbins();
    const int Nbins_3 = histo_etau_3->GetXaxis()->GetNbins();
    const int Nbins_4 = histo_etau_4->GetXaxis()->GetNbins();
    const int Nbins_5 = histo_etau_5->GetXaxis()->GetNbins();
    const int Nbins_6 = histo_etau_6->GetXaxis()->GetNbins();
  */
  if(histname == "Events_level_"){
    TH1F* hEvents_level = (TH1F*)((TH1F*)file_input->Get("Events_level"))->Clone(TString(name+"_"+histname));
    ofstream myfile;
    ofstream fs;
    if(sample_name=="data_obs"){
      myfile.open ("Event_number_mutau.csv");
      myfile << "Sample name" << "," << "Number used for scaling" << "," <<  "xsec"<< "," << "luminosity" << "," << "net weight"<< "," << "initial number of events"<< ","<< "events passing metFilter"<< "," << "events passing single ele trg" << "," << "Good ele"<< "," << "Good tau"<< ","<< "opposite charge selection"<< "," <<  "delta R cut"<< ","<<"third lepton veto"<< ","<<"bJet veto"<<","<< "Higgs pt cut" << "," <<"Visible mass cut"<<","<<  "events after met cut"<< std::endl;
    }
    else myfile.open ("Event_number_mutau.csv", ios::app);
        
    std::string SampleName = sample_name;
    std::string NumberEvents = std::to_string(ngen_);
    std::string CrossSection = std::to_string(xsec);
    std::string lumi = std::to_string(Luminosity);
    std::string NetWeight = std::to_string(net_weight);
    std::string FinalEvents_1 = std::to_string(hEvents_level->GetBinContent(1));
    std::string FinalEvents_2 = std::to_string(hEvents_level->GetBinContent(2));
    std::string FinalEvents_3 = std::to_string(hEvents_level->GetBinContent(3));
    std::string FinalEvents_4 = std::to_string(hEvents_level->GetBinContent(4));
    std::string FinalEvents_5 = std::to_string(hEvents_level->GetBinContent(5));
    std::string FinalEvents_6 = std::to_string(hEvents_level->GetBinContent(6));
    std::string FinalEvents_7 = std::to_string(hEvents_level->GetBinContent(7));
    std::string FinalEvents_8 = std::to_string(hEvents_level->GetBinContent(8));
    std::string FinalEvents_9 = std::to_string(hEvents_level->GetBinContent(9));
    std::string FinalEvents_10= std::to_string(hEvents_level->GetBinContent(10));
    std::string FinalEvents_11= std::to_string(hEvents_level->GetBinContent(11));
    std::string FinalEvents_12 = std::to_string(hEvents_level->GetBinContent(12));
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    
    myfile <<SampleName << ","  <<NumberEvents << "," << CrossSection << "," << lumi << "," << NetWeight <<","<<FinalEvents_1 <<","<< FinalEvents_2 <<","<<FinalEvents_3 << ","<<FinalEvents_4 <<","<<FinalEvents_5 <<","<<FinalEvents_6  <<","<<FinalEvents_7 <<","<< FinalEvents_8 <<","<<FinalEvents_9<<","<<FinalEvents_10 << ","<<FinalEvents_11 <<","<<FinalEvents_12  << std::endl;

    myfile.close();
    
  }

  outputFile->Close();

}


int main(int argc, char** argv)
{

  //std::cout<<"XXXXXXXXXXXXX This worked P1  XXXXXXXXXXXXX "<<'\n';
  
  std::string input = *(argv + 1);
  std::string output = *(argv + 2);
  std::string sample = *(argv + 3);
  std::string name = *(argv + 4);
  
  //std::string histNumber_1 = (*(argv + 6));
  //std::string histNumber_2 = (*(argv + 7));
  //int histNumber = stoi(HistNumber);
  //  cout<<"XXXXXXXXXXXXXXXXX   This works P2 histNumber_1 histNumber_2  "<< histNumber_1 <<"  " <<histNumber_2<<endl;
  float tes=0;
  if (argc > 1) {
    tes = atof(argv[5]);
  }
  TFile *f_Double = new TFile(input.c_str());
  cout<<"************* "<<" Input file         "<<"************ "<<endl;
  cout<<"XXXXXXXXXXXXX "<<input.c_str()<<"XXXXXXXXXXXXX "<<endl;
  cout<<"************* "<<" Output written to  "<<"************ "<<endl;
  cout<<"XXXXXXXXXXXXX "<<output.c_str()<<"XXXXXXXXXXXXX "<<endl;
  cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
  
  TH1F* nbevt = (TH1F*) f_Double->Get("Events_level" );
  //TTree *arbre = (TTree*) f_Double->Get("/ggNtuplizer/EventTree");
  //float ngen = 10000;
  std::cout<<"XXXXXXXXXXXXX This worked P1  XXXXXXXXXXXXX "<<'\n';

  //if (sample=="data_obs"){ Double_t ngen = nbevt->GetBinContent(2);}
  //else { Double_t ngen = nbevt->GetBinContent(1);}

  Double_t ngen = nbevt->GetBinContent(1);
  //  Double_t nFinal = nbevt->GetBinContent(12);

  std::cout<<"XXXXXXXXXXXXX  "<<"ngen = "<< ngen<<"   XXXXXXXXXXXXX "<<'\n';
  if (ngen<1)  std::cout<<"XXXXXXXXXXXXX  "<<" check this one !!" << input.c_str()<<"  XXXXXXXXXXXXX "<<'\n';

  Double_t xs=1.0; Double_t weight=1.0; Double_t luminosity=Luminosity;
  Double_t LOtoNNLO_DY = 5765.4/4954.0;
  Double_t LOtoNNLO_Wjets = 61526.7/50380; 
  if (sample=="ZL" || sample=="ZTT" || sample=="ZJ" || sample=="ZLL"){ xs=LOtoNNLO_DY*4954.0; weight=luminosity*xs/ngen;}
  else if (sample=="DY_LO" || sample=="DY_LOext"){ xs=5765.4; weight=luminosity*xs/ngen;}
  else if (sample=="ZTT1"){ xs=LOtoNNLO_DY*1012.5; weight=luminosity*xs/ngen;}
  else if (sample=="ZTT2"){ xs=LOtoNNLO_DY*332.8; weight=luminosity*xs/ngen;}
  else if (sample=="ZTT3"){ xs=LOtoNNLO_DY*101.8; weight=luminosity*xs/ngen;}
  else if (sample=="ZTT4"){ xs=LOtoNNLO_DY*54.8; weight=luminosity*xs/ngen;}
  else if (sample=="TTJets") {xs=831.76; weight=luminosity*xs/ngen;}
  else if (sample=="TTTo2L2Nu" ) {xs=88.29; weight=luminosity*xs/ngen;}
  else if (sample=="TTToHadronic" ) {xs=377.96; weight=luminosity*xs/ngen;}
  else if (sample=="TTToSemiLeptonic" ) {xs=365.35; weight=luminosity*xs/ngen;}

  /*  else if (sample=="ZTTinc"){ xs=LOtoNNLO_DY*5765.4; weight=2.458677925;}
  else if (sample=="ZTT1jet"){ xs=LOtoNNLO_DY*1012.5; weight=0.893891808;}
  else if (sample=="ZTT2jet"){ xs=LOtoNNLO_DY*332.8; weight=0.990764241;}
  else if (sample=="ZTT3jet"){ xs=LOtoNNLO_DY*101.8; weight=1.562190668;}
  else if (sample=="ZTT4jet"){ xs=LOtoNNLO_DY*54.8; weight=0.488914871;}
  */
  else if (sample=="ZTTinc"){ xs=LOtoNNLO_DY*4954.0; weight=4.877;}
  else if (sample=="ZTT1jet"){ xs=LOtoNNLO_DY*1012.5; weight=1.09;}
  else if (sample=="ZTT2jet"){ xs=LOtoNNLO_DY*332.8; weight=1.23;}
  else if (sample=="ZTT3jet"){ xs=LOtoNNLO_DY*101.8; weight=2.28;}
  else if (sample=="ZTT4jet"){ xs=LOtoNNLO_DY*54.8; weight=0.54;}
  else if (sample=="WJetsToLNu_inc") {xs=61526.7; weight=57.29;}
  else if (sample=="W1JetsToLNu") {xs=11778.3638; weight=7.81;}
  else if (sample=="W2JetsToLNu") {xs=3849.21974; weight=17.047;}
  else if (sample=="W3JetsToLNu") {xs=1166.04787; weight=2.36;}
  else if (sample=="W4JetsToLNu") {xs=593.055246; weight=2.098;}

  //else if (sample=="WJetsToLNu_inc") {xs=61526.7; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_2J" || sample=="WJetsToLNu"){xs=LOtoNNLO_Wjets*50380; weight=luminosity*xs/ngen;} // {xs=50380.0; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT100To200") {xs=LOtoNNLO_Wjets*1345.0; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT200To400") {xs=LOtoNNLO_Wjets*359.7; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT400To600") {xs=LOtoNNLO_Wjets*48.91; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT600To800") {xs=LOtoNNLO_Wjets*12.04; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT800To1200") {xs=LOtoNNLO_Wjets*5.52; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT1200To2500") {xs=LOtoNNLO_Wjets*1.33; weight=luminosity*xs/ngen;}
  else if (sample=="WJetsToLNu_HT2500ToInf") {xs=LOtoNNLO_Wjets*0.0322; weight=luminosity*xs/ngen;}

  else if (sample=="GluGluHToTauTau") {xs=44.14*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="VBFHToTauTau") {xs=3.782*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="ZHToTauTau") {xs=0.884*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="WplusH125" || sample=="WPlusHToTauTau") {xs=0.052; weight=luminosity*xs/ngen;}
  else if (sample=="WminusH125" || sample=="WMinusHToTauTau") {xs=0.0334; weight=luminosity*xs/ngen;}
  else if (sample=="GluGluHToWWTo2L2Nu") {xs=1.001; weight=luminosity*xs/ngen;}
  else if (sample=="VBFHToWWTo2L2Nu") {xs=0.0858; weight=luminosity*xs/ngen;}

  else if (sample=="QCD") {xs=720648000*0.00042; weight=luminosity*xs/ngen;}
  else if (sample=="data_obs"){weight=1.0;}


  else if (sample=="dataset_1"){weight=1.0;}
  else if (sample=="dataset_2"){weight=1.0;}
  else if (sample=="dataset_3"){weight=1.0;}
  else if (sample=="dataset_4"){weight=1.0;}
  else if (sample=="dataset_5"){weight=1.0;} 
   
 
  else if (sample=="ZJetsToNuNu_HT100To200") {xs=280.92; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT200To400") {xs=77.64; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT400To600") {xs=10.671; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT600To800") {xs=2.5611; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT800To1200") {xs=1.1778; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT1200To2500") {xs=0.2874; weight=luminosity*xs/ngen;}
  else if (sample=="ZJetsToNuNu_HT2500ToInf") {xs=0.006933; weight=luminosity*xs/ngen;}
  
  else if (sample=="WZTo1L3Nu") {xs=3.05; weight=luminosity*xs/ngen;}
  else if (sample=="WZTo1L1Nu2Q") {xs=10.71; weight=luminosity*xs/ngen;}
  else if (sample=="WZTo2L2Q") {xs=5.595; weight=luminosity*xs/ngen;}
  else if (sample=="WZTo3LNu") {xs=4.43; weight=luminosity*xs/ngen;}
  
  else if (sample=="ST_tW_antitop") {xs=35.6; weight=luminosity*xs/ngen;}
  else if (sample=="ST_tW_top") {xs=35.6; weight=luminosity*xs/ngen;}
  else if (sample=="ST_t-channel_antitop") {xs=80.95*3*0.105; weight=luminosity*xs/ngen;}
  else if (sample=="ST_t-channel_top") {xs=136.02*3*0.105; weight=luminosity*xs/ngen;}
  
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

  else if (sample=="ggH125") {xs=48.58*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="VBF125") {xs=3.782*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="ggH120") {xs=52.22*0.0698; weight=luminosity*xs/ngen;}
  else if (sample=="VBF120") {xs=3.935*0.0698; weight=luminosity*xs/ngen;}
  else if (sample=="ggH130") {xs=45.31*0.0541; weight=luminosity*xs/ngen;}
  else if (sample=="VBF130") {xs=3.637*0.0541; weight=luminosity*xs/ngen;}
  else if (sample=="ggH110") {xs=57.90*0.0791; weight=luminosity*xs/ngen;}
  else if (sample=="VBF110") {xs=4.434*0.0791; weight=luminosity*xs/ngen;}
  else if (sample=="ggH140") {xs=36.0*0.0360; weight=luminosity*xs/ngen;}
  else if (sample=="VBF140") {xs=3.492*0.0360; weight=luminosity*xs/ngen;}
  else if (sample=="ggH_WW125") {xs=48.58*0.2137*0.3258; weight=luminosity*xs/ngen;}
  else if (sample=="VBF_WW125") {xs=3.782*0.2137*0.3258; weight=luminosity*xs/ngen;}
  else if (sample=="WplusH120") {xs=0.9558*0.0698; weight=luminosity*xs/ngen;}
  else if (sample=="WplusH130") {xs=0.7414*0.0541; weight=luminosity*xs/ngen;}
  else if (sample=="WplusH110") {xs=1.335*0.0791; weight=luminosity*xs/ngen;}
  else if (sample=="WplusH140") {xs=0.6308*0.0360; weight=luminosity*xs/ngen;}
  else if (sample=="WminusH120") {xs=0.6092*0.0698; weight=luminosity*xs/ngen;}
  else if (sample=="WminusH130") {xs=0.4676*0.0541; weight=luminosity*xs/ngen;}
  else if (sample=="WminusH110") {xs=0.8587*0.0791; weight=luminosity*xs/ngen;}
  else if (sample=="WminusH140") {xs=0.394*0.0360; weight=luminosity*xs/ngen;}
  else if (sample=="ZH120") {xs=0.9939*0.0698; weight=luminosity*xs/ngen;}
  else if (sample=="ZH125") {xs=0.8839*0.0627; weight=luminosity*xs/ngen;}
  else if (sample=="ZH130") {xs=0.7899*0.0541; weight=luminosity*xs/ngen;}
  else if (sample=="ZH110") {xs=1.309*0.0791; weight=luminosity*xs/ngen;}
  else if (sample=="ZH140") {xs=0.6514*0.0360; weight=luminosity*xs/ngen;}
  else if (sample=="WGLNu") {xs=489.0; weight=luminosity*xs/ngen;}
  else if (sample=="WGstarMuMu") {xs=2.793; weight=luminosity*xs/ngen;}
  else if (sample=="WGstarEE") {xs=3.526; weight=luminosity*xs/ngen;}

  else if (sample=="EWKWMinus" || sample=="EWKWMinus2Jets") {xs=20.25; weight=luminosity*xs/ngen;}
  else if (sample=="EWKWPlus" || sample=="EWKWPlus2Jets") {xs=25.62; weight=luminosity*xs/ngen;}
  else if (sample=="EWKZLL" || sample=="EWKZ2Jets_ZToLL" || sample=="EWKZLL_TT" || sample=="EWKZLL_J" || sample=="EWKZLL_L" || sample=="EWKZLL_LL") {xs=3.987; weight=luminosity*xs/ngen;}
  else if (sample=="EWKZNuNu" || sample=="EWKZ2Jets_ZToNuNu" || sample=="EWKZNuNu_TT" || sample=="EWKZNuNu_J" || sample=="EWKZNuNu_L" || sample=="EWKZNuNu_LL") {xs=10.01; weight=luminosity*xs/ngen;}

  else {
    cout<<"Attention!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    cout<<"***********                                  *****************"<<endl;
    cout<<"***********"<<sample.c_str()<<"  NOT found in sample list ************"<<endl;
    cout<<"***********                                  *****************"<<endl;
    cout<<"Attention!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  }
  
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(10);
  
  
  std::vector<string> histnames;
  histnames.clear();
  std::vector<Double_t> leg_xoffsets;
  leg_xoffsets.clear();
  std::vector<Double_t> leg_yoffsets;
  leg_yoffsets.clear();
  std::vector<TString> xaxis_titles;
  xaxis_titles.clear();
  std::vector<TString> plotnames;
  plotnames.clear();

//  histnames.push_back(TString("_4"));
//  leg_xoffsets.push_back(0.);
//  leg_yoffsets.push_back(0.);
//  xaxis_titles.push_back(TString(""));
//  plotnames.push_back(TString(""));

  histnames.push_back("Muon_Pt_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("elePt"));

  histnames.push_back("Muon_En_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon #it{energy} [GeV]"));
  plotnames.push_back(TString("eleEn"));

  histnames.push_back("Muon_SCeta_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon SC  #eta"));
  plotnames.push_back(TString("eleSCEta"));

  histnames.push_back("Muon_SCphi_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon SC #phi"));
  plotnames.push_back(TString("eleSCPhi"));

  histnames.push_back("Muon_eta_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon #eta"));
  plotnames.push_back(TString("eleEta"));

  histnames.push_back("Muon_phi_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon #phi"));
  plotnames.push_back(TString("elePhi"));

  histnames.push_back("Tau_En_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Tau #it{energy} [GeV]"));
  plotnames.push_back(TString("tauEn"));

  histnames.push_back("Tau_Pt_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Tau #it{Pt} [GeV]"));
  plotnames.push_back(TString("tauPt"));

  histnames.push_back("Tau_eta_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Tau #eta"));
  plotnames.push_back(TString("tauEta"));

  histnames.push_back("Tau_phi_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Tau #phi"));
  plotnames.push_back(TString("tauPhi"));

  histnames.push_back("Tau_iso_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Tau isolation"));
  plotnames.push_back(TString("tauIso"));

  histnames.push_back("Tau_mass_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Tau mass"));
  plotnames.push_back(TString("tauMass"));

  histnames.push_back("Tau_Decay_Mode_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Tau Decay Mode"));
  plotnames.push_back(TString("tauDecayMode"));

  histnames.push_back("pfMET_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("pfMET [GeV]"));
  plotnames.push_back(TString("pfMET"));

  histnames.push_back("nJet_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Number of Jets"));
  plotnames.push_back(TString("nJet"));

  histnames.push_back("h_dPhi_");
  leg_xoffsets.push_back(-0.2);
  leg_yoffsets.push_back(-0.1);
  xaxis_titles.push_back(TString("#Delta#phi(electron, tau)"));
  plotnames.push_back(TString("dPhiEleTau"));

  histnames.push_back("Mt_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon-MET M_{T} [GeV]"));
  plotnames.push_back(TString("eleMETmT"));

  histnames.push_back("VisibleMass_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("VisibleMass_ [GeV]"));
  plotnames.push_back(TString("VisibleMass"));

  histnames.push_back("HiggsPt_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("HiggsPt_ [GeV]"));
  plotnames.push_back(TString("HiggsPt"));


  histnames.push_back("nVtx_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("nVtx_ [GeV]"));
  plotnames.push_back(TString("nVtx"));

  histnames.push_back("leadingJetPt_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("leadingJetPt_ [GeV]"));
  plotnames.push_back(TString("leadingJetPt"));

  histnames.push_back("Events_level_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Events_level_"));
  plotnames.push_back(TString("Evenets_level_"));

  histnames.push_back("MET_0_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("MET_0_"));
  plotnames.push_back(TString("MET_0_"));

  histnames.push_back("Cutflow_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Cutflow_"));
  plotnames.push_back(TString("Cutflow_"));

  cout<<"************************ works till here, after label *********************"<<endl;

  /* histnames.push_back("MET_1_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("MET_1_"));
  plotnames.push_back(TString("MET_1_"));

  histnames.push_back("MET_2_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("MET_2_"));
  plotnames.push_back(TString("MET_2_"));

  histnames.push_back("MET_3_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("MET_3_"));
  plotnames.push_back(TString("MET_3_"));

  histnames.push_back("MET_4_");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("MET_4_"));
  plotnames.push_back(TString("MET_4_"));
*/


  // cout<<"This works too P2" << endl;

  /*
  for(int i = 0; i < histnames.size(); i++){
    plot(histnames[i],leg_xoffsets[i],leg_yoffsets[i],xaxis_titles[i],plotnames[i]);
  }
  */
  TString output_ = TString(output);

  TFile* f_output;
  f_output = new TFile(output_ , "RECREATE");
  TDirectory *ptEle = f_output->mkdir("Muon_Pt_");
  TDirectory *enEle = f_output->mkdir("Muon_En_");
  TDirectory *phiEle = f_output->mkdir("Muon_phi_");
  TDirectory *etaEle = f_output->mkdir("Muon_eta_");
  TDirectory *SCPhiEle = f_output->mkdir("Muon_SCphi_");
  TDirectory *SCEtaEle = f_output->mkdir("Muon_SCeta_");
  TDirectory *ptTau = f_output->mkdir("Tau_Pt_");
  TDirectory *enTau = f_output->mkdir("Tau_En_");
  TDirectory *phiTau = f_output->mkdir("Tau_phi_");
  TDirectory *etaTau = f_output->mkdir("Tau_eta_");

  TDirectory *isoTau = f_output->mkdir("Tau_iso_");
  TDirectory *massTau = f_output->mkdir("Tau_mass_");
  TDirectory *dmTau = f_output->mkdir("Tau_Decay_Mode_");

  TDirectory *pfMET = f_output->mkdir("pfMET_");
  TDirectory *h_dPhi = f_output->mkdir("h_dPhi_");
  TDirectory *Mt = f_output->mkdir("Mt_");
  TDirectory *VisibleMass = f_output->mkdir("VisibleMass_");
  TDirectory *HiggsPt = f_output->mkdir("HiggsPt_");
  TDirectory *nVtx = f_output->mkdir("nVtx_");
  TDirectory *leadingJetPt = f_output->mkdir("leadingJetPt_");
  TDirectory *nJet = f_output->mkdir("nJet_");
  TDirectory *Events_level_all = f_output->mkdir("Events_level_");
  TDirectory *Cutflow_all = f_output->mkdir("Cutflow_"); 
  TDirectory *MET_0_ = f_output->mkdir("MET_0_");
  // TDirectory *MET_1_ = f_output->mkdir("MET_1_");
  //TDirectory *MET_2_ = f_output->mkdir("MET_2_");
  //TDirectory *MET_3_ = f_output->mkdir("MET_3_");
  //TDirectory *MET_4_ = f_output->mkdir("MET_4_");


  f_output->Close();

  // cout<<"This workks too before makehist" <<endl;
  make_hist(input, output, "Muon_En_", sample, name, weight, true, ngen , xs);
  make_hist(input, output, "Tau_En_", sample , name,weight, true, ngen , xs);
  make_hist(input, output, "Muon_eta_", sample ,  name,weight, true, ngen , xs);
  make_hist(input, output, "Tau_eta_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "Muon_Pt_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "Tau_Pt_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "Muon_phi_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "Tau_phi_", sample ,name,weight,  true, ngen , xs);

  make_hist(input, output, "Tau_iso_", sample ,name,weight,  true, ngen , xs); 
  make_hist(input, output, "Tau_mass_", sample ,name,weight,  true, ngen , xs); 
  make_hist(input, output, "Tau_Decay_Mode_", sample ,name,weight,  true, ngen , xs); 

  make_hist(input, output, "nVtx_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "pfMET_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "h_dPhi_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "nJet_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "leadingJetPt_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "VisibleMass_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "HiggsPt_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "Mt_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "Events_level_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "Cutflow_", sample ,name,weight, true, ngen , xs);
  make_hist(input, output, "MET_0_", sample ,name,weight, true, ngen , xs);
  //make_hist(input, output, "MET_1_", sample ,name,weight, true, ngen , xs);
  //make_hist(input, output, "MET_2_", sample ,name,weight, true, ngen , xs);
  //make_hist(input, output, "MET_3_", sample ,name,weight, true, ngen , xs);
  //make_hist(input, output, "MET_4_", sample ,name,weight, true, ngen , xs);

  // cout<<"This works too, after makehist" <<endl;

}



