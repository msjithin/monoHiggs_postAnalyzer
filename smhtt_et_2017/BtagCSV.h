#ifndef BTAGCSV_h
#define BTAGCSV_h

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

struct BTagCSV {
  struct BTagSF {
    int op;
    TString measurement;
    TString sys;
    int jetFlavor;
    float etaMin,etaMax,ptMin,ptMax,discrMin,discrMax;
    TF1 formula;
    BTagSF(TString csvline);
    inline TString GetName() { return TString( std::to_string(op) )+"_"+measurement+"_"+sys+"_"+TString( std::to_string(jetFlavor) ); }
    float EvalSF(float pt, float eta);
  };
  std::map<TString,BTagSF*> sfmap;
  
  BTagCSV(TString csvname);
  BTagCSV::BTagSF* getBTagSF(int op=1, TString measurement="comb", TString sys="central", int jetFlavor=0);
  float EvalSF(int op, TString measurement, TString sys, int jetFlavor, float pt, float eta);
};

std::vector<TString> split(TString str,TString delim);

#endif
