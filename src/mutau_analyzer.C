/////mutau_analyzer.C
// For use with Ntuples made from ggNtuplizer
// Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
// To compile using rootcom to an executable named 'analyze':
//$ ./rootcom mutau_analyzer analyze
//
// To run, assuming this is compiled to an executable named 'analyze':
//$ ./analyze /hdfs/store/user/jmadhusu/LatestNtuples/ /afs/hep.wisc.edu/user/ms/CMSSW_9_4_4/src/2017_analysis/etau/output.root -1 10000
//./analyze /hdfs/store/user/jmadhusu/MonoHiggs_MC2017/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/crab_ZZZ/180603_185329/0000/ /afs/hep.wisc.edu/user/ms/CMSSW_9_4_4/src/2017_analysis/analyzer/output.root -1 10000
// Runs over every event in the folder LatestNtuples, reporting progress every 10000 events
// and storing the resulting histograms in the file output.root.
//
// To plot, for example, single photon trigger efficiency as a function of photon pt:
//$ root -l
// root[0] TFile *f = new TFile("output.root");
// root[1] TGraphAsymmErrors *efficiency = new TGraphAsymmErrors((TH1F*)f->Get("Photon_Et_300_2"),(TH1F*)f->Get("Photon_Et_300_1"));
// root[2] efficiency->Draw("AP")
// root[3] efficiency->SetTitle("Single photon trigger efficiency")
// root[4] efficiency->GetXaxis()->SetTitle("Photon p_{T}")
// root[5] efficiency->Draw("AP")
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
#include "object_functions.h"
#include "fill_histograms.h"
#include "selections.h"
#include "get_met_systenatics.h"

using namespace std;
using std::vector;
int main(int argc, const char *argv[])
{
  TStopwatch sw;
  sw.Start();

  myMap1 = new map<string, TH1F *>();
  myMap2 = new map<string, TH2F *>();
  std::string SampleName = argv[7];
  std::string isMC = argv[6];
  std::string outputfile = argv[2];
  Long64_t maxEvents = atof(argv[3]);
  string sp = "0";
  cout << "argc = " << argc << endl;
  if (argc > 8)
    sp = string(argv[8]);
  cout << "sp = " << sp << endl;
  const char *signalpara = sp.c_str();

  if (maxEvents < -1LL)
  {
    std::cout << "Please enter a valid value for maxEvents (parameter 3)." << std::endl;
    return 1;
  }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
  {
    std::cout << "Please enter a valid value for reportEvery (parameter 4) " << std::endl;
    return 1;
  }
  // std::string SampleName = argv[5];

  mutau_analyzer t(argv[1], argv[2], isMC, SampleName, signalpara);
  t.Loop(maxEvents, reportEvery, SampleName);
  // delete myMap1;
  cout << " Outpt written to " << outputfile << endl;
  sw.Stop();
  sw.Print();
  return 0;
}

void mutau_analyzer::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
{

  int nTotal;
  nTotal = 0;
  int report_ = 0;
  int report_test = 0;
  double numberOfEvents = 0;
  int nInspected;
  nInspected = 0;
  double nInspected_genWeighted;
  nInspected_genWeighted = 0.0;
  debug = false;
  if (fChain == 0)
    return;

  std::vector<int> muCand;
  muCand.clear();
  std::vector<int> tauCand;
  tauCand.clear();
  std::vector<int> aisrtauCand;
  aisrtauCand.clear();
  TString sample = TString(SampleName);

  int nHiggs = 0;
  int nHToMuTau = 0;
  int found_mt = 0;
  int muCand_1 = 0;
  int muCand_2 = 0;
  int muCand_3 = 0;
  int tauCand_1 = 0;
  int tauCand_2 = 0;
  int tauCand_3 = 0;

  Double_t Pt_Bins[26] = {0.0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
  Double_t Pt_Bins_highPt[21] = {100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};

  // TH1F* h_cutflow=new TH1F("cutflow", "cutflow", 10, 0, 10); h_cutflow->Sumw2();
  TH1F *h_cutflow_n = new TH1F("cutflow_n", "cutflow_n", 8, 0, 8);
  h_cutflow_n->Sumw2();
  TH1F *h_cutflow_n_fr = new TH1F("cutflow_n_fr", "cutflow_n_fr", 8, 0, 8);
  h_cutflow_n_fr->Sumw2();
  TH1F *h_cutflow_n_dyll = new TH1F("cutflow_n_dyll", "cutflow_n_dyll", 8, 0, 8);
  h_cutflow_n_dyll->Sumw2();
  TH1F *h_cutflow_n_dyll_fr = new TH1F("cutflow_n_dyll_fr", "cutflow_n_dyll_fr", 8, 0, 8);
  h_cutflow_n_dyll_fr->Sumw2();
  // TH1F* h_cutflow_Htt=new TH1F("cutflow_Htt", "cutflow_Htt", 11, 0, 11); h_cutflow_Htt->Sumw2();

  // if(debug)cout<<" setting up other files ..."<<endl;
  // RoccoR  rc("sf_files/roCorr_Run2_v3/RoccoR2017.txt");

  Long64_t nentries = fChain->GetEntries();
  if (is_MC == true)
    std::cout << ".... MC file ..... " << std::endl;
  else
    std::cout << ".... DATA file ..... " << std::endl;

  std::cout << "Coming in: " << std::endl;
  std::cout << "nentries:" << nentries << std::endl;
  // Look at up to maxEvents events, or all if maxEvents == -1.
  Long64_t nentriesToCheck = nentries;
  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;

  std::cout << "Running over " << nTotal << " events." << std::endl;
  // TStopwatch sw;
  // sw.Start();

  for (Long64_t jentry = 0; jentry < nentriesToCheck; jentry++)
  {
    if (debug)
      cout << "event " << jentry << endl;
    muCand.clear();
    tauCand.clear();
    aisrtauCand.clear();
    muCand_nom.clear();
    tauCand_nom.clear();
    aisrtauCand_nom.clear();
    jetCand.clear();
    jetCand_nom.clear();

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    double inspected_event_weight = 1.0;
    if (is_MC)
      fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight / fabs(genWeight) : inspected_event_weight = 0.0;
    nInspected_genWeighted += inspected_event_weight;
    nInspected += 1;
    // h_insEvents->SetBinContent(1, nInspected_genWeighted);
    //=1.0 for real data
    double event_weight = 1.0;
    double weight = 1.0;
    double pileup_sf = 1.0;
    double applySf = 1.0;
    bool passMetTrigger = false;
    bool passCrossTrigger = false;
    int report_i = 0;
    bool Ztt_selector = false;

    numberOfEvents += weight;
    if (is_MC)
      weight = inspected_event_weight;
    else
      weight = 1.0;
    if (is_MC)
      pileup_sf = h_pileup->GetBinContent(h_pileup->GetXaxis()->FindBin(puTrue->at(0)));
    weight = weight * pileup_sf;
    if (is_MC && prefiringweight!=0)
      weight = weight * prefiringweight;
    if (isGoodVtx == false)
      continue;
    t_index = get_t_Cand();
    tbar_index = get_tbar_Cand();

    /////Trigger bit selection
    if ( HLTMet >> 0 & 1 == 1 || HLTMet >> 1 & 1 == 1)
      passMetTrigger = true;

    ////

    /////
    if (debug)
      cout << "entry # : " << jentry << endl;

    if (!is_MC)
      event_weight = 1.0;
    else
      event_weight = weight;

    // //cout<< "zprimeBaryonic_signal "<<zprimeBaryonic_signal<<"  "<<__LINE__<<endl;
    // if (found_ZprimeBaryonic==true && std::find(signalParameters->begin(), signalParameters->end(), zprimeBaryonic_signal) == signalParameters->end()) {
    // 	//cout<< "zprimeBaryonic_signal "<<zprimeBaryonic_signal<<"  "<<endl;
    // 	continue;
    // }
    // if (found_ZprimeBaryonic==true )
    // 	plotFill("nEvents_ZpB", zprimeBaryonic_signal, 50, 0, 50,  1.0);

    if (metFilters == 0 && (passMetTrigger))
    {
      eventNumber = jentry;
      nSingleTrgPassed += event_weight;
      // if(check_unc)cout<<"Preparing NOMINAL"<<endl;
      orginal_jetPt.clear();
      for (float pt : (*jetPt))
        orginal_jetPt.push_back(pt);

      selections(event_weight, 0, "nominal"); // this is for nominal

      // // jet-tau fakes
      // selections(event_weight,  1, "jetFakes");
      // selections(event_weight,  -1, "jetFakes");
      // if(is_MC){
      //   /// UP
      //   string shape_names[7] = {
      //     //"tauES",
      //     //"JES","JER",
      //     //"metresponse", "metresolution",  "metunclustered",
      //     "tauIDunc", "tauTRGunc", "leptonTRGunc",
      //     "prefiringUnc", "muonMissID",
      //     "dyShape", "ttbarShape"};
      //   for (int i = 0; i < 7 ; i++ ){
      //     // if (shape_names[i] != "metunclustered")
      //     // 	continue;
      //     // cout<<"shape name = "<<shape_names[i]<<endl;
      //     selections(event_weight,  1, shape_names[i]);
      //     selections(event_weight,  -1, shape_names[i]);
      //   }
      // }
    }

    report_test = nentriesToCheck / 20;
    while (report_test > 10)
    {
      report_test = report_test / 10;
      report_i++;
    }
    if (nentriesToCheck > 20)
      reportEvery = report_test * pow(10, report_i);
    else
      reportEvery = 1;
    if (jentry % reportEvery == 0)
    {
      std::cout << "Finished entry " << jentry << "/" << (nentriesToCheck - 1) << std::endl;
    }
  }

  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  if ((nentriesToCheck - 1) % reportEvery != 0)
    std::cout << "Finished entry " << (nentriesToCheck - 1) << "/" << (nentriesToCheck - 1) << std::endl;
  // sw.Stop();
  std::cout << "All events checked." << std::endl;
  std::cout << "*******************************************" << std::endl;
  std::cout << "******************Jithin's original*************************" << std::endl;
  std::cout << std::setw(20) << std::right << "Initial entries " << numberOfEvents << std::endl;
  std::cout << std::setw(20) << std::right << "Passing smikking " << nPassedSkimmed << std::endl;
  std::cout << std::setw(20) << std::right << "Inspected genWeightd " << nInspected_genWeighted << std::setw(10) << std::right << "   % change= " << (numberOfEvents - nInspected_genWeighted) * 100 / numberOfEvents << std::endl;
  std::cout << std::setw(20) << std::right << "GoodMuonPassed " << nGoodMuonPassed << std::setw(10) << std::right << "   % change= " << (nSingleTrgPassed - nGoodMuonPassed) * 100 / nSingleTrgPassed << std::endl;
  std::cout << std::setw(20) << std::right << "GoodTauPassed " << nGoodTauPassed << std::setw(10) << std::right << "   % change= " << (nGoodMuonPassed - nGoodTauPassed) * 100 / nGoodMuonPassed << std::endl;
  //   std::cout<<std::setw(20) <<std::right <<"TauIsoPassed "<<nTauIsoPassed<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nTauIsoPassed)*100/nGoodTauPassed<<std::endl;
  // std::cout<<std::setw(20) <<std::right <<"TauDecayModePassed "<<nTauDecayModePassed<<std::setw(10) <<std::right << "   % change= "<<(nTauIsoPassed-nTauDecayModePassed)*100/nTauIsoPassed<<std::endl;

  std::cout << std::setw(20) << std::right << "opp charge " << nGoodMuTauPassed << std::setw(10) << std::right << "   % change= " << (nGoodTauPassed - nGoodMuTauPassed) * 100 / nGoodTauPassed << std::endl;

  std::cout << std::setw(20) << std::right << "PassedThirdLepVeto " << nPassedThirdLepVeto << std::setw(10) << std::right << "   % change= " << (nGoodMuTauPassed - nPassedThirdLepVeto) * 100 / nGoodMuTauPassed << std::endl;
  std::cout << std::setw(20) << std::right << "PassedBjetVeto " << nPassedBjetVeto << std::setw(10) << std::right << "   % change= " << (nPassedThirdLepVeto - nPassedBjetVeto) * 100 / nPassedThirdLepVeto << std::endl;
  std::cout << std::setw(20) << std::right << "DeltaRPassed " << nDeltaRPassed << std::setw(10) << std::right << "   % change= " << (nPassedBjetVeto - nDeltaRPassed) * 100 / nPassedBjetVeto << std::endl;

  std::cout << std::setw(20) << std::right << "Total change :" << (numberOfEvents - nDeltaRPassed) * 100 / numberOfEvents << std::endl;
  std::cout << "*******************************************" << std::endl;
  std::cout << "*******************************************" << std::endl;
  std::cout << std::setw(20) << std::right << "Number of events inspected: " << nInspected << std::endl;
  std::cout << std::setw(20) << std::right << "Number of events inspected (minus negative gen. weights): " << nInspected_genWeighted << std::endl;

  h_cutflow_n->SetBinContent(1, nInspected_genWeighted);
  h_cutflow_n->SetBinContent(2, nSingleTrgPassed);
  h_cutflow_n->SetBinContent(3, nGoodMuonPassed);
  h_cutflow_n->SetBinContent(4, nGoodTauPassed);
  h_cutflow_n->SetBinContent(5, nGoodMuTauPassed);
  h_cutflow_n->SetBinContent(6, nPassedThirdLepVeto);
  h_cutflow_n->SetBinContent(7, nPassedBjetVeto);
  h_cutflow_n->SetBinContent(8, nDeltaRPassed);

  h_cutflow_n_fr->SetBinContent(1, nInspected_genWeighted);
  h_cutflow_n_fr->SetBinContent(2, nSingleTrgPassed_fr);
  h_cutflow_n_fr->SetBinContent(3, nGoodMuonPassed_fr);
  h_cutflow_n_fr->SetBinContent(4, nGoodTauPassed_fr);
  h_cutflow_n_fr->SetBinContent(5, nGoodMuTauPassed_fr);
  h_cutflow_n_fr->SetBinContent(6, nPassedThirdLepVeto_fr);
  h_cutflow_n_fr->SetBinContent(7, nPassedBjetVeto_fr);
  h_cutflow_n_fr->SetBinContent(8, nDeltaRPassed_fr);

  /// dy Z->ll
  h_cutflow_n_dyll->SetBinContent(1, nInspected_genWeighted);
  h_cutflow_n_dyll->SetBinContent(2, nSingleTrgPassed_dyll);
  h_cutflow_n_dyll->SetBinContent(3, nGoodMuonPassed_dyll);
  h_cutflow_n_dyll->SetBinContent(4, nGoodTauPassed_dyll);
  h_cutflow_n_dyll->SetBinContent(5, nGoodMuTauPassed_dyll);
  h_cutflow_n_dyll->SetBinContent(6, nPassedThirdLepVeto_dyll);
  h_cutflow_n_dyll->SetBinContent(7, nPassedBjetVeto_dyll);
  h_cutflow_n_dyll->SetBinContent(8, nDeltaRPassed_dyll);

  h_cutflow_n_dyll_fr->SetBinContent(1, nInspected_genWeighted);
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
  cout<<"------------------"<<endl;
  cout<<"nboosted"<<"   "<<"nresolved"<<endl;
  cout<<nboosted<<"          "<<nresolved<<endl;
  cout<<"                  "<<endl;
}

void mutau_analyzer::BookHistos(const char *file1, const char *file2)
{
  TFile *file_in = new TFile(file1, "READ");
  fileName = new TFile(file2, "RECREATE");

  // makeOutputTree(tree);
  fileName->cd();
  // cout<<"cloning nEvents hist"<<endl;
  h_nEvents = (TH1F *)((TH1F *)file_in->Get("nEvents"))->Clone(TString("nEvents"));
  file_in->Close();
}

// Fill the sequential histos at a particular spot in the sequence

void mutau_analyzer::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{

  //  h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);
}

void mutau_analyzer::printP4values(string when = "")
{
  if (check_unc == false)
    return;
  printTabSeparated("entry # : ", eventNumber, "before/after = ", when, "\n",
                    "Shift ", unc_shift, "\n",
                    "electron pt", my_muP4.Pt(),
                    "electron eta", my_muP4.Eta(), "\n",
                    "tau pt", my_tauP4.Pt(),
                    "tau eta", my_tauP4.Eta(), "\n",
                    "MET", my_metP4.Pt(), "\n");
}
