#include "HTTutilities/Jet2TauFakes/interface/FakeFactor.h"

void test(TString fname="/afs/cern.ch/user/m/mspanrin/public/Htautau/FakeRate/et/incl/fakeFactors_20180412_tight.root"){

  // Retrieve the fake factor
  TFile* ff_file = TFile::Open(fname);
  FakeFactor* ff    = (FakeFactor*)ff_file->Get("ff_comb");

  // Fill inputs
  std::vector<string> inputNames( ff->inputs() ) ;
  std::vector<double> inputs( inputNames.size() );
  inputs[0] = 30; //tau_pt;
  inputs[1] = 0;  //tau_decayMode;
  inputs[2] = 1;  //njet
  inputs[3] = 40; //mvis;
  inputs[4] = 10; //mt;
  inputs[5] = 0.00; //muon_iso;
  inputs[6] = 0.4;  //frac_qcd
  inputs[7] = 0.3;  //frac_w
  inputs[8] = 0.2;  //frac_tt

  // Retrieve fake factors
  double ff_nom = ff->value(inputs); // nominal fake factor

  double syst_qcd_up = ff->value(inputs, "ff_qcd_syst_up");
  double syst_qcd_down = ff->value(inputs, "ff_qcd_syst_down");
  double syst_w_up = ff->value(inputs, "ff_w_syst_up");
  double syst_w_down = ff->value(inputs, "ff_w_syst_down");
  double syst_tt_up = ff->value(inputs, "ff_tt_syst_up");
  double syst_tt_down = ff->value(inputs, "ff_tt_syst_down");

  double stat_qcd_up = ff->value(inputs, "ff_qcd_stat_up");
  double stat_qcd_down = ff->value(inputs, "ff_qcd_stat_down");
  double stat_w_up = ff->value(inputs, "ff_w_stat_up");
  double stat_w_down = ff->value(inputs, "ff_w_stat_down");
  double stat_tt_up = ff->value(inputs, "ff_tt_stat_up");
  double stat_tt_down = ff->value(inputs, "ff_tt_stat_down");
  
  for( int i = 0; i < inputNames.size(); i++ ){
    cout << inputNames[i]<< "= " << inputs[i] << ", ";
  }
  cout <<std::endl <<  "ff= " << ff_nom << endl;
  cout << " ----- Systematic uncertainties ----- " << endl;
  cout << "Uncertainties on corrections: " << endl;
  cout << "syst(tt)= " << (syst_tt_up-ff_nom)/ff_nom << "%, syst(w+dy)= " << (syst_w_up-ff_nom)/ff_nom*100 << "%, syst(qcd)= " << (syst_qcd_up-ff_nom)/ff_nom*100 << "%" <<  endl;
  cout << "Uncertainties on fake factors: " << endl;
  cout << "stat(tt)= " << (stat_tt_up-ff_nom)/ff_nom*100 << "%, stat(w+dy)= " << (stat_w_up-ff_nom)/ff_nom*100 << "%, stat(qcd)= " << (stat_qcd_up-ff_nom)/ff_nom*100 << "%" << endl;

  delete ff;
  ff_file->Close();

}
