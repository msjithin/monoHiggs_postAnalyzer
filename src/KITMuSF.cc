#ifndef KIT_MU_SF_h
#define KIT_MU_SF_h

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TString.h>
#include <TRandom.h>
#include <TMath.h>
#include <assert.h>
#include "string.h"


using namespace std;
class KITMuSF {
private:
  TFile* fileIn;
  TH1* etaBinsH;
  int nEtaBins ;
  bool debug_;
  //map<string, TGraph*>* myMap_effData;
  //map<string, TGraph*>* myMap_effMC;
  //int findPtBin(double pt);
  double get_EfficiencyData(double pt, double eta);
  double get_EfficiencyMC(double pt, double eta);
  TGraph* ZMassEtaLt0p9_Data;
  TGraph* ZMassEta0p9to1p2_Data;
  TGraph* ZMassEta1p2to2p1_Data;
  TGraph* ZMassEtaGt2p1_Data;
  TGraph* ZMassEtaLt0p9_MC;
  TGraph* ZMassEta0p9to1p2_MC;
  TGraph* ZMassEta1p2to2p1_MC;
  TGraph* ZMassEtaGt2p1_MC;
public:
  KITMuSF();
  ~KITMuSF();
  void init_ScaleFactors( string fileName);
  double get_ScaleFactor(double pt, double eta);
};
KITMuSF::KITMuSF(){
  debug_=false;
  if(debug_)cout<<"KITMuSF constr"<<endl;  
  
}
KITMuSF::~KITMuSF(){
  if(debug_)cout<<"KITMuSF destr"<<endl;
  fileIn->Close();
}
void KITMuSF::init_ScaleFactors(string fileName){
  //if(debug_)cout<<"inputRootFile : "<<inputRootFile<<endl;
  TFile * file = new TFile(fileName.c_str());
  if (file->IsZombie()) {
    std::cout << "file " << fileName << " is not found...   quitting " << std::endl;
    exit(-1);
  }
  fileIn = TFile::Open(fileName.c_str());
  etaBinsH=(TH1*) fileIn->Get("etaBinsH");
  etaBinsH->SetDirectory(0);
  nEtaBins = etaBinsH->GetNbinsX();
  ZMassEtaLt0p9_Data = (TGraph*) fileIn->Get("ZMassEtaLt0p9_Data");
  ZMassEta0p9to1p2_Data = (TGraph*) fileIn->Get("ZMassEta0p9to1p2_Data");
  ZMassEta1p2to2p1_Data = (TGraph*) fileIn->Get("ZMassEta1p2to2p1_Data");
  ZMassEtaGt2p1_Data = (TGraph*) fileIn->Get("ZMassEtaGt2p1_Data");
  ZMassEtaLt0p9_MC = (TGraph*) fileIn->Get("ZMassEtaLt0p9_MC");
  ZMassEta0p9to1p2_MC = (TGraph*) fileIn->Get("ZMassEta0p9to1p2_MC");
  ZMassEta1p2to2p1_MC = (TGraph*) fileIn->Get("ZMassEta1p2to2p1_MC");
  ZMassEtaGt2p1_MC = (TGraph*) fileIn->Get("ZMassEtaGt2p1_MC");
  
}
double KITMuSF::get_EfficiencyData(double pt, double eta){
  double rv_eff= 1.0;
  if(pt >100) pt=99;
  else if (pt<=10) pt=10;
  if     (abs(eta)<0.9) rv_eff = ZMassEtaLt0p9_Data->Eval(pt);
  else if(abs(eta)<1.2) rv_eff = ZMassEta0p9to1p2_Data->Eval(pt);
  else if(abs(eta)<2.1) rv_eff = ZMassEta1p2to2p1_Data->Eval(pt);
  else                  rv_eff = ZMassEtaGt2p1_Data->Eval(pt);
  if(debug_)cout<<"pt "<<pt<<" eta "<<eta<<" eff_d "<<rv_eff<<endl; 
  return rv_eff;
}
double KITMuSF::get_EfficiencyMC(double pt, double eta){
  double rv_eff= 1.0;
  if(pt >100) pt=99;
  else if (pt<=10) pt=10;
  if     (abs(eta)<0.9) rv_eff = ZMassEtaLt0p9_MC->Eval(pt);
  else if(abs(eta)<1.2) rv_eff = ZMassEta0p9to1p2_MC->Eval(pt);
  else if(abs(eta)<2.1) rv_eff = ZMassEta1p2to2p1_MC->Eval(pt);
  else                  rv_eff = ZMassEtaGt2p1_MC->Eval(pt);
  if(debug_)cout<<"pt "<<pt<<" eta "<<eta<<" eff_mc "<<rv_eff<<endl;
  return rv_eff;
}
double KITMuSF::get_ScaleFactor(double pt, double eta){
  double rv_sf = 1.0;
  double eff_data = get_EfficiencyData(pt, eta);
  double eff_mc = get_EfficiencyMC(pt, eta);
  if (eff_mc !=0)
    rv_sf = eff_data / eff_mc ;
  else
    rv_sf = 1.0;
  if(debug_)cout<<"pt "<<pt<<" eta "<<eta<<" rv_sf "<<rv_sf<<endl;
  return rv_sf;
}
#endif
