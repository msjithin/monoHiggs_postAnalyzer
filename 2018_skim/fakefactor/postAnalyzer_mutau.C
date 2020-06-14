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
//typedef std::vector<double> NumV;

using namespace std;
//Double_t Luminosity = 41521.0;//Lumi for inclusive
//Double_t Luminosity = 4823.0;    // for run B
//Double_t Luminosity = 9664.0;  // for run C
//Double_t Luminosity = 4252.0;  // for run D
//Double_t Luminosity = 9278.0;  // for run E
//Double_t Luminosity = 13540.0; // for run F



int main(int argc, char** argv)
{

  //std::cout<<"XXXXXXXXXXXXX This worked P1  XXXXXXXXXXXXX "<<'\n';
  bool debug = false;
  std::string input_file = *(argv + 1);
  std::string output_file = *(argv + 2);
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
  //TFile *f_Double = new TFile(input.c_str());
  cout<<"************* "<<" Input file         "<<" ************ "<<endl;
  cout<<"************* "<<input_file.c_str()         <<" ************ "<<endl;
  cout<<"************* "<<" Output written to  "<<" ************ "<<endl;
  cout<<"*************  "<<output_file.c_str()       <<" ************ "<<endl;
  cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
  
  if (debug==true) std::cout<<"XXXXXXXXXXXXX This worked P1  XXXXXXXXXXXXX "<<'\n';
  TFile *file_input = new TFile(TString(input_file));
  TString output_ = TString(output_file);

  TFile *f_output = new TFile(output_ , "RECREATE");
  if (debug==true) std::cout<<"XXXXXXXXXXXXX This worked P2  XXXXXXXXXXXXX "<<'\n';
  

  //  TH1F* histo_tauPt_test  = (TH1F*)(file_input->Get("tauPt_inc_test"))->Clone("tauPt_inc_test");
  TH1F* histo_tauPt_t_dm0  = (TH1F*)(file_input->Get("tauPt_dm0_tight"))->Clone("tauPt_dm0_tight");
  if (debug==true) std::cout<<"XXXXXXXXXXXXX This worked P2.a  XXXXXXXXXXXXX "<<'\n';

  TH1F* histo_tauPt_t_dm1  = (TH1F*)(file_input->Get("tauPt_dm1_tight"))->Clone("tauPt_dm1_tight");
  if (debug==true) std::cout<<"XXXXXXXXXXXXX This worked P2.b  XXXXXXXXXXXXX "<<'\n';

  TH1F* histo_tauPt_t_dm10  = (TH1F*)(file_input->Get("tauPt_dm_10_tight"))->Clone("tauPt_dm_10_tight");
  TH1F* histo_tauPt_t_dm11  = (TH1F*)(file_input->Get("tauPt_dm_11_tight"))->Clone("tauPt_dm_11_tight");

  if (debug==true) std::cout<<"XXXXXXXXXXXXX This worked P2.c  XXXXXXXXXXXXX "<<'\n';

  TH1F* histo_tauPt_vl_dm0 = (TH1F*)(file_input->Get("tauPt_dm0_veryloose"))->Clone("tauPt_dm0_veryloose");
  TH1F* histo_tauPt_vl_dm1 = (TH1F*)(file_input->Get("tauPt_dm1_veryloose"))->Clone("tauPt_dm1_veryloose");
  TH1F* histo_tauPt_vl_dm10 = (TH1F*)(file_input->Get("tauPt_dm10_veryloose"))->Clone("tauPt_dm10_veryloose");
  TH1F* histo_tauPt_vl_dm11 = (TH1F*)(file_input->Get("tauPt_dm11_veryloose"))->Clone("tauPt_dm11_veryloose");

  // TH1F* histo_tauPt_m_dm0 = (TH1F*)(file_input->Get("tauPt_dm0_medium"))->Clone("tauPt_dm0_medium");
  // TH1F* histo_tauPt_m_dm1 = (TH1F*)(file_input->Get("tauPt_dm1_medium"))->Clone("tauPt_dm1_medium");
  // TH1F* histo_tauPt_m_dm10 = (TH1F*)(file_input->Get("tauPt_dm_10_medium"))->Clone("tauPt_dm_10_medium");

  Float_t pt_bins[8]= {20, 25, 30, 35, 40, 50, 60, 120};
  if (debug==true) std::cout<<"XXXXXXXXXXXXX This worked P4  XXXXXXXXXXXXX "<<'\n';

  TH1F* h_tauFF_dm0= new TH1F("tauFF_dm0","tau fake factor dm0",7,pt_bins ); h_tauFF_dm0->Sumw2();
  TH1F* h_tauFF_dm1= new TH1F("tauFF_dm1","tau fake factor dm1",7,pt_bins ); h_tauFF_dm1->Sumw2();
  TH1F* h_tauFF_dm10= new TH1F("tauFF_dm10","tau fake factor dm10",7,pt_bins ); h_tauFF_dm10->Sumw2();
  

  h_tauFF_dm0->Add(histo_tauPt_t_dm0);
  h_tauFF_dm0->Divide(histo_tauPt_vl_dm0);
  h_tauFF_dm1->Add(histo_tauPt_t_dm1);
  h_tauFF_dm1->Divide(histo_tauPt_vl_dm1);
  h_tauFF_dm10->Add(histo_tauPt_t_dm10);
  h_tauFF_dm10->Divide(histo_tauPt_vl_dm10);

  h_tauFF_dm0->SetMinimum(0.0);
  if (h_tauFF_dm0->GetMaximum() > 0.5) 
    h_tauFF_dm0->SetMaximum(h_tauFF_dm0->GetMaximum()*1.10);
  else
    h_tauFF_dm0->SetMaximum(0.5);

  h_tauFF_dm1->SetMinimum(0.0);
  if (h_tauFF_dm1->GetMaximum() > 0.5) 
    h_tauFF_dm1->SetMaximum(h_tauFF_dm1->GetMaximum()*1.10);
  else
    h_tauFF_dm1->SetMaximum(0.5);

  h_tauFF_dm10->SetMinimum(0.0);
  if (h_tauFF_dm10->GetMaximum() > 0.5) 
    h_tauFF_dm10->SetMaximum(h_tauFF_dm10->GetMaximum()*1.10);
  else
    h_tauFF_dm10->SetMaximum(0.5);

  
  TGraphAsymmErrors* g_tauFF_dm0 = new TGraphAsymmErrors();
  g_tauFF_dm0->Divide(histo_tauPt_t_dm0, histo_tauPt_vl_dm0);
  g_tauFF_dm0->GetXaxis()->SetTitle("tau pt[GeV]");
  g_tauFF_dm0->GetYaxis()->SetTitle("jet-tau fake rate");
  g_tauFF_dm0->SetMinimum(0.0);
  g_tauFF_dm0->SetMaximum(0.5);
  g_tauFF_dm0->SetTitle("tau fake factor decay mode 0");
  g_tauFF_dm0->Write("hpt_dm0_tight_hpt_dm0_veryloose");
  
  TGraphAsymmErrors* g_tauFF_dm1 = new TGraphAsymmErrors();
  g_tauFF_dm1->Divide(histo_tauPt_t_dm1, histo_tauPt_vl_dm1);
  g_tauFF_dm1->GetXaxis()->SetTitle("tau pt[GeV]");
  g_tauFF_dm1->GetYaxis()->SetTitle("jet-tau fake rate");
  g_tauFF_dm1->SetMinimum(0.0);
  g_tauFF_dm1->SetMaximum(0.5);
  g_tauFF_dm1->SetTitle("tau fake factor decay mode 1");
  g_tauFF_dm1->Write("hpt_dm1_tight_hpt_dm1_veryloose");

  TGraphAsymmErrors* g_tauFF_dm10 = new TGraphAsymmErrors();
  g_tauFF_dm10->Divide(histo_tauPt_t_dm10, histo_tauPt_vl_dm10);
  g_tauFF_dm10->GetXaxis()->SetTitle("tau pt[GeV]");
  g_tauFF_dm10->GetYaxis()->SetTitle("jet-tau fake rate");
  g_tauFF_dm10->SetMinimum(0.0);
  g_tauFF_dm10->SetMaximum(0.5);
  g_tauFF_dm10->SetTitle("tau fake factor decay mode 10");
  g_tauFF_dm10->Write("hpt_dm10_tight_hpt_dm10_veryloose");

  TGraphAsymmErrors* g_tauFF_dm11 = new TGraphAsymmErrors();
  g_tauFF_dm11->Divide(histo_tauPt_t_dm11, histo_tauPt_vl_dm11);
  g_tauFF_dm11->GetXaxis()->SetTitle("tau pt[GeV]");
  g_tauFF_dm11->GetYaxis()->SetTitle("jet-tau fake rate");
  g_tauFF_dm11->SetMinimum(0.0);
  g_tauFF_dm11->SetMaximum(0.5);
  g_tauFF_dm11->SetTitle("tau fake factor decay mode 11");
  g_tauFF_dm11->Write("hpt_dm11_tight_hpt_dm11_veryloose");

  // TGraphAsymmErrors* g_tauFF_m_dm0 = new TGraphAsymmErrors();
  // g_tauFF_m_dm0->Divide(histo_tauPt_m_dm0, histo_tauPt_vl_dm0);
  // g_tauFF_m_dm0->GetXaxis()->SetTitle("tau pt[GeV]");
  // g_tauFF_m_dm0->GetYaxis()->SetTitle("jet-tau fake rate");
  // g_tauFF_m_dm0->SetMinimum(0.0);
  // g_tauFF_m_dm0->SetMaximum(0.5);
  // g_tauFF_m_dm0->SetTitle("tau fake factor decay mode 0");
  // g_tauFF_m_dm0->Write("hpt_dm0_medium_hpt_dm0_veryloose");
  // TGraphAsymmErrors* g_tauFF_m_dm1 = new TGraphAsymmErrors();
  // g_tauFF_m_dm1->Divide(histo_tauPt_m_dm1, histo_tauPt_vl_dm1);
  // g_tauFF_m_dm1->GetXaxis()->SetTitle("tau pt[GeV]");
  // g_tauFF_m_dm1->GetYaxis()->SetTitle("jet-tau fake rate");
  // g_tauFF_m_dm1->SetMinimum(0.0);
  // g_tauFF_m_dm1->SetMaximum(0.5);
  // g_tauFF_m_dm1->SetTitle("tau fake factor decay mode 1");
  // g_tauFF_m_dm1->Write("hpt_dm1_medium_hpt_dm1_veryloose");
  // TGraphAsymmErrors* g_tauFF_m_dm10 = new TGraphAsymmErrors();
  // g_tauFF_m_dm10->Divide(histo_tauPt_m_dm10, histo_tauPt_vl_dm10);
  // g_tauFF_m_dm10->GetXaxis()->SetTitle("tau pt[GeV]");
  // g_tauFF_m_dm10->GetYaxis()->SetTitle("jet-tau fake rate");
  // g_tauFF_m_dm10->SetMinimum(0.0);
  // g_tauFF_m_dm10->SetMaximum(0.5);
  // g_tauFF_m_dm10->SetTitle("tau fake factor decay mode 10");
  // g_tauFF_m_dm10->Write("hpt_dm10_medium_hpt_dm10_veryloose");

  //  histo_tauPt_test->Write();
  h_tauFF_dm0->Write();
  h_tauFF_dm1->Write();
  h_tauFF_dm10->Write();
  histo_tauPt_t_dm0->Write();
  histo_tauPt_t_dm1->Write();
  histo_tauPt_t_dm10->Write();
  histo_tauPt_vl_dm0->Write();
    
  if (debug==true) std::cout<<"XXXXXXXXXXXXX This worked P7  XXXXXXXXXXXXX "<<'\n';
  file_input->Close();
  f_output->Close();

}



