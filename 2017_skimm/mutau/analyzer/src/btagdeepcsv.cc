#include "btagdeepcsv.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>


btagdeepcsv::btagdeepcsv(const std::string& str) : fName(str)
{
  ifstream csvfile(fName);
  if ( !csvfile.is_open() ) {
    cout << "Unable to read " << fName << endl;
    return;
  }
  else {
    cout << "File read " << fName << endl;
    std::string csvline;
    std::getline(csvfile,csvline);
    while ( std::getline(csvfile,csvline) ) {
      BTagSF* btagsf = new BTagSF(csvline);
      sfmap[btagsf->GetName()] = btagsf;
      
    }
  }
}
btagdeepcsv::~btagdeepcsv()
{
}
btagdeepcsv::BTagSF(TString csvline) {
  vector<TString> values = split(csvline,",");
  op = values[0].Atoi();
  measurement = values[1].ReplaceAll(" ","");
  sys = values[2].ReplaceAll(" ","");
  jetFlavor = values[3].Atoi();
  etaMin = values[4].Atof();
  etaMax = values[5].Atof();
  ptMin = values[6].Atof();
  ptMax = values[7].Atof();
  discrMin = values[8].Atoi();
  discrMax = values[9].Atoi();
  formula = TF1(GetName(),values[10].ReplaceAll(" ","").ReplaceAll("\"",""),ptMin,ptMax);
}

float BTagCSV::BTagSF::EvalSF(float pt,float eta) {
  if ( fabs(eta) < etaMin || fabs(eta) > etaMax ) return 1;
  return formula.Eval(pt);
}

BTagCSV::BTagCSV(TString csvname) {
  ifstream csvfile(csvname);
  if ( !csvfile.is_open() ) {
    cout << "Unable to read " << csvname << endl;
    return;
  }
  std::string csvline;
  std::getline(csvfile,csvline);
  while ( std::getline(csvfile,csvline) ) {
    BTagSF* btagsf = new BTagSF(csvline);
    sfmap[btagsf->GetName()] = btagsf;
  }
}

BTagCSV::BTagSF* BTagCSV::getBTagSF(int op, TString measurement, TString sys, int jetFlavor) {
  TString sfname = TString( std::to_string(op) )+"_"+measurement+"_"+sys+"_"+TString( std::to_string(jetFlavor) );
  return sfmap[sfname];
}

float BTagCSV::EvalSF(int op, TString measurement, TString sys, int jetFlavor, float pt, float eta) {
  return getBTagSF(op,measurement,sys,jetFlavor)->EvalSF(pt,eta);
}

vector<TString> split(TString str,TString delim) {
  vector<TString> splitTString;
  TString token;
  int idx = 0;
  while ( str.Tokenize(token,idx,delim) ) splitTString.push_back(token);
  return splitTString;
}
