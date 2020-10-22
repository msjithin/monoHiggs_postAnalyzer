#ifndef BTAGDEEPCSV_h
#define BTAGDEEPCSV_h

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TString.h>
#include <TRandom.h>
#include <TMath.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <map>


class btagdeepcsv {
  std::string fName;
 public:
  btagdeepcsv(const std::string& str);
  ~btagdeepcsv();
  
  int op;
  TString measurement;
  TString sys;
  int jetFlavor;
  float etaMin,etaMax,ptMin,ptMax,discrMin,discrMax;
  TF1 formula;
  void BTagSF(TString csvline);
  inline TString GetName() { return TString( std::to_string(op) )+"_"+measurement+"_"+sys+"_"+TString( std::to_string(jetFlavor) ); }
  float EvalSF(float pt, float eta);
  
  std::map<TString,BTagSF*> sfmap;
  
  void BTagCSV(TString csvname);
  double getBTagSF(int op=1, TString measurement="comb", TString sys="central", int jetFlavor=0);
  float EvalSF(int op, TString measurement, TString sys, int jetFlavor, float pt, float eta);

  std::vector<TString> split(TString str,TString delim);
  
 
};


#endif
