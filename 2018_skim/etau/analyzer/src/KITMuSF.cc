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
  TGraph* ZMassEtaLt1p0_Data;
  TGraph* ZMassEta1p0to1p48_Data;
  TGraph* ZMassEta1p48to1p65_Data;
  TGraph* ZMassEta1p65to2p1_Data;
  TGraph* ZMassEtaGt2p1_Data;
  TGraph* ZMassEtaLt1p0_MC;
  TGraph* ZMassEta1p0to1p48_MC;
  TGraph* ZMassEta1p48to1p65_MC;
  TGraph* ZMassEta1p65to2p1_MC;
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
  ZMassEtaLt1p0_Data = (TGraph*) fileIn->Get("ZMassEtaLt1p0_Data");
  ZMassEta1p0to1p48_Data = (TGraph*) fileIn->Get("ZMassEta1p0to1p48_Data");
  ZMassEta1p48to1p65_Data = (TGraph*) fileIn->Get("ZMassEta1p48to1p65_Data");
  ZMassEta1p65to2p1_Data = (TGraph*) fileIn->Get("ZMassEta1p65to2p1_Data");
  ZMassEtaGt2p1_Data = (TGraph*) fileIn->Get("ZMassEtaGt2p1_Data");
  ZMassEtaLt1p0_MC = (TGraph*) fileIn->Get("ZMassEtaLt1p0_MC");
  ZMassEta1p0to1p48_MC = (TGraph*) fileIn->Get("ZMassEta1p0to1p48_MC");
  ZMassEta1p48to1p65_MC = (TGraph*) fileIn->Get("ZMassEta1p48to1p65_MC");
  ZMassEta1p65to2p1_MC = (TGraph*) fileIn->Get("ZMassEta1p65to2p1_MC");
  ZMassEtaGt2p1_MC = (TGraph*) fileIn->Get("ZMassEtaGt2p1_MC");

}
double KITMuSF::get_EfficiencyData(double pt, double eta){
  double rv_eff= 1.0;
  if(pt >1000) pt=999;
  else if (pt<=10) pt=10;
  if     (abs(eta)<1.0) rv_eff = ZMassEtaLt1p0_Data->Eval(pt);
  else if(abs(eta)<1.48) rv_eff = ZMassEta1p0to1p48_Data->Eval(pt);
  else if(abs(eta)<1.65) rv_eff = ZMassEta1p48to1p65_Data->Eval(pt);
  else if(abs(eta)<2.1) rv_eff = ZMassEta1p65to2p1_Data->Eval(pt);
  else                  rv_eff = ZMassEtaGt2p1_Data->Eval(pt);
  if(debug_)cout<<"pt "<<pt<<" eta "<<eta<<" eff_d "<<rv_eff<<endl; 
  return rv_eff;
}
double KITMuSF::get_EfficiencyMC(double pt, double eta){
  double rv_eff= 1.0;
  if(pt >1000) pt=999;
  else if (pt<=10) pt=10;
  if     (abs(eta)<1.0) rv_eff = ZMassEtaLt1p0_MC->Eval(pt);
  else if(abs(eta)<1.48) rv_eff = ZMassEta1p0to1p48_MC->Eval(pt);
  else if(abs(eta)<1.65) rv_eff = ZMassEta1p48to1p65_MC->Eval(pt);
  else if(abs(eta)<2.1) rv_eff = ZMassEta1p65to2p1_MC->Eval(pt);
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
