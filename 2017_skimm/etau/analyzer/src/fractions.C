#ifndef FRANCTIONS_C
#define FRANCTIONS_C




#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"

//#include "myevent.h"
//#include "LinkDef.h"
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
using namespace std;



void getFractions(int category, double mvis, double& frac_qcd, double& frac_w, double& frac_tt ){
  
  if( category==1 ) // 0 jet , ptH<10
    {
      if(mvis<50){frac_qcd=0; frac_w=0; frac_tt=0;}
      if(mvis>50 && mvis < 100)  {frac_qcd=0.86; frac_w=0.14; frac_tt=0;}
      if(mvis>100 && mvis < 150) {frac_qcd=0.76; frac_w=0.24; frac_tt=0;}
      if(mvis>150 ) {frac_qcd=0.61; frac_w=0.39; frac_tt=0;}
    }

  else if( category==2 ) // 0 jet , ptH>10
    {
      if(mvis<50){frac_qcd=0.43; frac_w=0.57; frac_tt=0;}
      if(mvis>50 && mvis < 100)  {frac_qcd=0.5; frac_w=0.5; frac_tt=0;}
      if(mvis>100 && mvis < 150) {frac_qcd=0.49; frac_w=0.51; frac_tt=0;}
      if(mvis>150 && mvis < 200) {frac_qcd=0.39; frac_w=0.61; frac_tt=0;}
      if(mvis>200 ) {frac_qcd=0.38; frac_w=0.61; frac_tt=0;}
    }
  else if( category==3 ) // boosted 1 jet
    {
      if(mvis<50){frac_qcd=0.68; frac_w=0.31; frac_tt=0.01;}
      if(mvis>50 && mvis < 100)  {frac_qcd=0.52; frac_w=0.47; frac_tt=0.01;}
      if(mvis>100 && mvis < 150) {frac_qcd=0.48; frac_w=0.51; frac_tt=0.01;}
      if(mvis>150 && mvis < 200) {frac_qcd=0.46; frac_w=0.53; frac_tt=0.01;}
      if(mvis>200 ) {frac_qcd=0.48; frac_w=0.50; frac_tt=0.02;}
    }
  else if( category==4 ) // boosted >=2 jet
    {
      if(mvis<50){frac_qcd=0.52; frac_w=0.44; frac_tt=0.04;}
      if(mvis>50 && mvis < 100)  {frac_qcd=0.24; frac_w=0.68; frac_tt=0.08;}
      if(mvis>100 && mvis < 150) {frac_qcd=0.20; frac_w=0.70; frac_tt=0.10;}
      if(mvis>150 && mvis < 200) {frac_qcd=0.14; frac_w=0.74; frac_tt=0.12;}
      if(mvis>200 ) {frac_qcd=0.13; frac_w=0.67; frac_tt=0.10;}
    }
  else if( category==5 ) // vbf ptH < 200
    {
      if(mvis<50){frac_qcd=0.50; frac_w=0.46; frac_tt=0.04;}
      if(mvis>50 && mvis < 100)  {frac_qcd=0.45; frac_w=0.50; frac_tt=0.05;}
      if(mvis>100 && mvis < 150) {frac_qcd=0.34; frac_w=0.60; frac_tt=0.06;}
      if(mvis>150 && mvis < 200) {frac_qcd=0.22; frac_w=0.69; frac_tt=0.09;}
      if(mvis>200 ) {frac_qcd=0.26; frac_w=0.66; frac_tt=0.08;}
    } 
  else if( category==6 ) // vbf ptH > 200
    {
      if(mvis<50)                {frac_qcd=0.73; frac_w=0.01; frac_tt=0.26;}
      if(mvis>50 && mvis < 100)  {frac_qcd=0.0; frac_w=0.81; frac_tt=0.19;}
      if(mvis>100 && mvis < 150) {frac_qcd=0.0; frac_w=0.88; frac_tt=0.12;}
      if(mvis>150 && mvis < 200) {frac_qcd=0.0; frac_w=0.82; frac_tt=0.18;}
      if(mvis>200 ) {frac_qcd=0.0; frac_w=0.82; frac_tt=0.18;}
    }

}


#endif  //FRANCTIONS_C
