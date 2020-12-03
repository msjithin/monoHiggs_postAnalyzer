#include "TF1.h"
#include "TFile.h"
#include "TString.h"

using namespace std;

class applyFF_w_lpt
{
public:
  TFile* frawff;
  TF1* ff_qcd_0jet;;
  TF1* ff_qcd_1jet;
  TF1* ff_qcd_2jet;
  TF1* ff_w_0jet;
  TF1* ff_w_1jet;
  TF1* ff_w_2jet;
  TF1* ff_tt_0jet;

  TF1* ff_qcd_0jet_up1;
  TF1* ff_qcd_1jet_up1;
  TF1* ff_qcd_2jet_up1;
  TF1* ff_w_0jet_up1;
  TF1* ff_w_1jet_up1;
  TF1* ff_w_2jet_up1;
  TF1* ff_tt_0jet_up1;
  TF1* ff_qcd_0jet_up2;
  TF1* ff_qcd_1jet_up2;
  TF1* ff_qcd_2jet_up2;
  TF1* ff_w_0jet_up2;
  TF1* ff_w_1jet_up2;
  TF1* ff_w_2jet_up2;
  TF1* ff_tt_0jet_up2;
  TF1* ff_qcd_0jet_down1;
  TF1* ff_qcd_1jet_down1;
  TF1* ff_qcd_2jet_down1;
  TF1* ff_w_0jet_down1;
  TF1* ff_w_1jet_down1;
  TF1* ff_w_2jet_down1;
  TF1* ff_tt_0jet_down1;
  TF1* ff_qcd_0jet_down2;
  TF1* ff_qcd_1jet_down2;
  TF1* ff_qcd_2jet_down2;
  TF1* ff_w_0jet_down2;
  TF1* ff_w_1jet_down2;
  TF1* ff_w_2jet_down2;
  TF1* ff_tt_0jet_down2;

  TFile* fmvisclosure; 
  TF1* lptclosure_w_taupt30to40;
  TF1* lptclosure_w_taupt40to50;
  TF1* lptclosure_w_tauptgt50;
  TF1* lptclosure_qcd_taupt30to40;
  TF1* lptclosure_qcd_taupt40to50;
  TF1* lptclosure_qcd_tauptgt50;
  TF1* lptclosure_tt_taupt30to40;
  TF1* lptclosure_tt_taupt40to50;
  TF1* lptclosure_tt_tauptgt50;
  TF1* lptclosure_w_xtrg_taupt30to40;
  TF1* lptclosure_w_xtrg_taupt40to50;
  TF1* lptclosure_w_xtrg_tauptgt50;
  TF1* lptclosure_qcd_xtrg_taupt30to40;
  TF1* lptclosure_qcd_xtrg_taupt40to50;
  TF1* lptclosure_qcd_xtrg_tauptgt50;
  TF1* lptclosure_tt_xtrg_taupt30to40;
  TF1* lptclosure_tt_xtrg_taupt40to50;
  TF1* lptclosure_tt_xtrg_tauptgt50;


  TFile* fosssclosure ;
  TF1* osssclosure_qcd;
  TF1* mtclosure_w;
  
  TF1* osssclosure_qcd_up1;
  TF1* mtclosure_w_up1;
  TF1* osssclosure_qcd_up2;
  TF1* mtclosure_w_up2;
  TF1* osssclosure_qcd_down1;
  TF1* mtclosure_w_down1;
  TF1* osssclosure_qcd_down2;
  TF1* mtclosure_w_down2;
  
  applyFF_w_lpt(){
    //std::cout<<"applyFF_w_lpt constructed"<<endl;
    frawff = TFile::Open("sf_files/FF_lptCorrInBinsOfTauPt/ff_files_et_2018/uncorrected_fakefactors_et.root");
    ff_qcd_0jet=(TF1*) frawff->Get("rawFF_et_qcd_0jet");
    ff_qcd_1jet=(TF1*) frawff->Get("rawFF_et_qcd_1jet");
    ff_qcd_2jet=(TF1*) frawff->Get("rawFF_et_qcd_2jet");
    ff_w_0jet=(TF1*) frawff->Get("rawFF_et_w_0jet");
    ff_w_1jet=(TF1*) frawff->Get("rawFF_et_w_1jet");
    ff_w_2jet=(TF1*) frawff->Get("rawFF_et_w_2jet");
    ff_tt_0jet=(TF1*) frawff->Get("mc_rawFF_et_tt");

    ff_qcd_0jet_up1=(TF1*) frawff->Get("rawFF_et_qcd_0jet_unc1_up");
    ff_qcd_1jet_up1=(TF1*) frawff->Get("rawFF_et_qcd_1jet_unc1_up");
    ff_qcd_2jet_up1=(TF1*) frawff->Get("rawFF_et_qcd_2jet_unc1_up");
    ff_w_0jet_up1=(TF1*) frawff->Get("rawFF_et_w_0jet_unc1_up");
    ff_w_1jet_up1=(TF1*) frawff->Get("rawFF_et_w_1jet_unc1_up");
    ff_w_2jet_up1=(TF1*) frawff->Get("rawFF_et_w_2jet_unc1_up");
    ff_tt_0jet_up1=(TF1*) frawff->Get("mc_rawFF_et_tt_unc1_up");
    ff_qcd_0jet_up2=(TF1*) frawff->Get("rawFF_et_qcd_0jet_unc2_up");
    ff_qcd_1jet_up2=(TF1*) frawff->Get("rawFF_et_qcd_1jet_unc2_up");
    ff_qcd_2jet_up2=(TF1*) frawff->Get("rawFF_et_qcd_2jet_unc2_up");
    ff_w_0jet_up2=(TF1*) frawff->Get("rawFF_et_w_0jet_unc2_up");
    ff_w_1jet_up2=(TF1*) frawff->Get("rawFF_et_w_1jet_unc2_up");
    ff_w_2jet_up2=(TF1*) frawff->Get("rawFF_et_w_2jet_unc2_up");
    ff_tt_0jet_up2=(TF1*) frawff->Get("mc_rawFF_et_tt_unc2_up");
    ff_qcd_0jet_down1=(TF1*) frawff->Get("rawFF_et_qcd_0jet_unc1_down");
    ff_qcd_1jet_down1=(TF1*) frawff->Get("rawFF_et_qcd_1jet_unc1_down");
    ff_qcd_2jet_down1=(TF1*) frawff->Get("rawFF_et_qcd_2jet_unc1_down");
    ff_w_0jet_down1=(TF1*) frawff->Get("rawFF_et_w_0jet_unc1_down");
    ff_w_1jet_down1=(TF1*) frawff->Get("rawFF_et_w_1jet_unc1_down");
    ff_w_2jet_down1=(TF1*) frawff->Get("rawFF_et_w_2jet_unc1_down");
    ff_tt_0jet_down1=(TF1*) frawff->Get("mc_rawFF_et_tt_unc1_down");
    ff_qcd_0jet_down2=(TF1*) frawff->Get("rawFF_et_qcd_0jet_unc2_down");
    ff_qcd_1jet_down2=(TF1*) frawff->Get("rawFF_et_qcd_1jet_unc2_down");
    ff_qcd_2jet_down2=(TF1*) frawff->Get("rawFF_et_qcd_2jet_unc2_down");
    ff_w_0jet_down2=(TF1*) frawff->Get("rawFF_et_w_0jet_unc2_down");
    ff_w_1jet_down2=(TF1*) frawff->Get("rawFF_et_w_1jet_unc2_down");
    ff_w_2jet_down2=(TF1*) frawff->Get("rawFF_et_w_2jet_unc2_down");
    ff_tt_0jet_down2=(TF1*) frawff->Get("mc_rawFF_et_tt_unc2_down");

    fmvisclosure = TFile::Open("sf_files/FF_lptCorrInBinsOfTauPt/ff_files_et_2018/FF_corrections_1.root");
    lptclosure_w_taupt30to40=(TF1*) fmvisclosure->Get("closure_lpt_taupt30to40_et_w");
    lptclosure_w_taupt40to50=(TF1*) fmvisclosure->Get("closure_lpt_taupt40to50_et_w");
    lptclosure_w_tauptgt50=(TF1*) fmvisclosure->Get("closure_lpt_tauptgt50_et_w");
    lptclosure_qcd_taupt30to40=(TF1*) fmvisclosure->Get("closure_lpt_taupt30to40_et_qcd");
    lptclosure_qcd_taupt40to50=(TF1*) fmvisclosure->Get("closure_lpt_taupt40to50_et_qcd");
    lptclosure_qcd_tauptgt50=(TF1*) fmvisclosure->Get("closure_lpt_tauptgt50_et_qcd");
    lptclosure_tt_taupt30to40=(TF1*) fmvisclosure->Get("closure_lpt_taupt30to40_et_ttmc");
    lptclosure_tt_taupt40to50=(TF1*) fmvisclosure->Get("closure_lpt_taupt40to50_et_ttmc");
    lptclosure_tt_tauptgt50=(TF1*) fmvisclosure->Get("closure_lpt_tauptgt50_et_ttmc");
    lptclosure_w_xtrg_taupt30to40=(TF1*) fmvisclosure->Get("closure_lpt_taupt30to40_xtrg_et_w");
    lptclosure_w_xtrg_taupt40to50=(TF1*) fmvisclosure->Get("closure_lpt_taupt40to50_xtrg_et_w");
    lptclosure_w_xtrg_tauptgt50=(TF1*) fmvisclosure->Get("closure_lpt_tauptgt50_xtrg_et_w");
    lptclosure_qcd_xtrg_taupt30to40=(TF1*) fmvisclosure->Get("closure_lpt_taupt30to40_xtrg_et_qcd");
    lptclosure_qcd_xtrg_taupt40to50=(TF1*) fmvisclosure->Get("closure_lpt_taupt40to50_xtrg_et_qcd");
    lptclosure_qcd_xtrg_tauptgt50=(TF1*) fmvisclosure->Get("closure_lpt_tauptgt50_xtrg_et_qcd");
    lptclosure_tt_xtrg_taupt30to40=(TF1*) fmvisclosure->Get("closure_lpt_taupt30to40_xtrg_et_ttmc");
    lptclosure_tt_xtrg_taupt40to50=(TF1*) fmvisclosure->Get("closure_lpt_taupt40to50_xtrg_et_ttmc");
    lptclosure_tt_xtrg_tauptgt50=(TF1*) fmvisclosure->Get("closure_lpt_tauptgt50_xtrg_et_ttmc");


    fosssclosure  = TFile::Open("sf_files/FF_lptCorrInBinsOfTauPt/ff_files_et_2018/FF_QCDcorrectionOSSS.root");
    osssclosure_qcd=(TF1*) fosssclosure->Get("closure_OSSS_dr_flat_et_qcd");
    mtclosure_w=(TF1*) fosssclosure->Get("closure_mt_et_w");

    osssclosure_qcd_up1=(TF1*) fosssclosure->Get("closure_OSSS_mvis_et_qcd_unc1_up");
    mtclosure_w_up1=(TF1*) fosssclosure->Get("closure_mt_et_w_unc1_up");
    osssclosure_qcd_up2=(TF1*) fosssclosure->Get("closure_OSSS_mvis_et_qcd_unc2_up");
    mtclosure_w_up2=(TF1*) fosssclosure->Get("closure_mt_et_w_unc2_up");
    osssclosure_qcd_down1=(TF1*) fosssclosure->Get("closure_OSSS_mvis_et_qcd_unc1_down");
    mtclosure_w_down1=(TF1*) fosssclosure->Get("closure_mt_et_w_unc1_down");
    osssclosure_qcd_down2=(TF1*) fosssclosure->Get("closure_OSSS_mvis_et_qcd_unc2_down");
    mtclosure_w_down2=(TF1*) fosssclosure->Get("closure_mt_et_w_unc2_down");

  }
  ~applyFF_w_lpt(){
    //std::cout<<"applyFF_w_lpt destructing..."<<endl;
    frawff->Close();
    fmvisclosure->Close();
    fosssclosure->Close();
  }
  float get_raw_FF(float pt, TF1* fct){
    float ff=1.0;
    ff=fct->Eval(pt);
    if (pt>100) ff=fct->Eval(100);
    return ff;
    //if (ff<0) return 0.0;
    //else if (ff>2) return 2.0;
    //else return ff;
  }

  float get_mvis_closure(float mvis, TF1* fct){
    float corr=1.0;
    corr=fct->Eval(mvis);
    if (mvis>250) corr=fct->Eval(250);
    return corr;
    //if (corr<0 or corr>2) return 1.0;
    //else return corr;
  }

  float get_lpt_closure(float lpt, TF1* fct){
    float corr=1.0;
    corr=fct->Eval(lpt);
    if (lpt>120) corr=fct->Eval(120);
    return corr;
  }

  float get_mt_closure(float mt, TF1* fct){
    float corr=1.0;
    corr=fct->Eval(mt);
    if (corr<0 or corr>2) return 1.0;
    else return corr;
  }

  float get_ff(float pt, float mt, float mvis, float msv, float lpt, float met, int njets, bool xtrg, float frac_tt, float frac_qcd, float frac_w, TString shift){
    float ff_qcd=1.0;
    float ff_w=0;
    float ff_tt=1.0;

    // Raw FF
    if (njets==0){
      ff_qcd=get_raw_FF(pt,ff_qcd_0jet);
      if (shift=="ff_qcd_up1" or shift=="ff_qcd_up1_0jet") ff_qcd=get_raw_FF(pt,ff_qcd_0jet_up1);
      else if (shift=="ff_qcd_up2" or shift=="ff_qcd_up2_0jet") ff_qcd=get_raw_FF(pt,ff_qcd_0jet_up2);
      else if (shift=="ff_qcd_down1" or shift=="ff_qcd_down1_0jet") ff_qcd=get_raw_FF(pt,ff_qcd_0jet_down1);
      else if (shift=="ff_qcd_down2" or shift=="ff_qcd_down2_0jet") ff_qcd=get_raw_FF(pt,ff_qcd_0jet_down2);
      ff_w=get_raw_FF(pt,ff_w_0jet);
      if (shift=="ff_w_up1" or shift=="ff_w_up1_0jet") ff_w=get_raw_FF(pt,ff_w_0jet_up1);
      else if (shift=="ff_w_up2" or shift=="ff_w_up2_0jet") ff_w=get_raw_FF(pt,ff_w_0jet_up2);
      else if (shift=="ff_w_down1" or shift=="ff_w_down1_0jet") ff_w=get_raw_FF(pt,ff_w_0jet_down1);
      else if (shift=="ff_w_down2" or shift=="ff_w_down2_0jet") ff_w=get_raw_FF(pt,ff_w_0jet_down2);
    }
    else if (njets==1){
      ff_qcd=get_raw_FF(pt,ff_qcd_1jet);
      if (shift=="ff_qcd_up1" or shift=="ff_qcd_up1_1jet") ff_qcd=get_raw_FF(pt,ff_qcd_1jet_up1);
      else if (shift=="ff_qcd_up2" or shift=="ff_qcd_up2_1jet") ff_qcd=get_raw_FF(pt,ff_qcd_1jet_up2);
      else if (shift=="ff_qcd_down1" or shift=="ff_qcd_down1_1jet") ff_qcd=get_raw_FF(pt,ff_qcd_1jet_down1);
      else if (shift=="ff_qcd_down2" or shift=="ff_qcd_down2_1jet") ff_qcd=get_raw_FF(pt,ff_qcd_1jet_down2);
      ff_w=get_raw_FF(pt,ff_w_1jet);
      if (shift=="ff_w_up1" or shift=="ff_w_up1_1jet") ff_w=get_raw_FF(pt,ff_w_1jet_up1);
      else if (shift=="ff_w_up2" or shift=="ff_w_up2_1jet") ff_w=get_raw_FF(pt,ff_w_1jet_up2);
      else if (shift=="ff_w_down1" or shift=="ff_w_down1_1jet") ff_w=get_raw_FF(pt,ff_w_1jet_down1);
      else if (shift=="ff_w_down2" or shift=="ff_w_down2_1jet") ff_w=get_raw_FF(pt,ff_w_1jet_down2);
    }
    else{
      ff_qcd=get_raw_FF(pt,ff_qcd_2jet);
      if (shift=="ff_qcd_up1" or shift=="ff_qcd_up1_2jet") ff_qcd=get_raw_FF(pt,ff_qcd_2jet_up1);
      else if (shift=="ff_qcd_up2" or shift=="ff_qcd_up2_2jet") ff_qcd=get_raw_FF(pt,ff_qcd_2jet_up2);
      else if (shift=="ff_qcd_down1" or shift=="ff_qcd_down1_2jet") ff_qcd=get_raw_FF(pt,ff_qcd_2jet_down1);
      else if (shift=="ff_qcd_down2" or shift=="ff_qcd_down2_2jet") ff_qcd=get_raw_FF(pt,ff_qcd_2jet_down2);
      ff_w=get_raw_FF(pt,ff_w_2jet);
      if (shift=="ff_w_up1" or shift=="ff_w_up1_2jet") ff_w=get_raw_FF(pt,ff_w_2jet_up1);
      else if (shift=="ff_w_up2" or shift=="ff_w_up2_2jet") ff_w=get_raw_FF(pt,ff_w_2jet_up2);
      else if (shift=="ff_w_down1" or shift=="ff_w_down1_2jet") ff_w=get_raw_FF(pt,ff_w_2jet_down1);
      else if (shift=="ff_w_down2" or shift=="ff_w_down2_2jet") ff_w=get_raw_FF(pt,ff_w_2jet_down2);
    }
    ff_tt=get_raw_FF(pt,ff_tt_0jet);
    if (shift=="ff_tt_up1") ff_tt=get_raw_FF(pt,ff_tt_0jet_up1);
    else if (shift=="ff_tt_up2") ff_tt=get_raw_FF(pt,ff_tt_0jet_up2);
    else if (shift=="ff_tt_down1") ff_tt=get_raw_FF(pt,ff_tt_0jet_down1);
    else if (shift=="ff_tt_down2") ff_tt=get_raw_FF(pt,ff_tt_0jet_down2);

    // L pT closure
    if (xtrg){
      if (pt<40){
	if (shift=="lptclosure_xtrg_qcd_up") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_xtrg_taupt30to40)*1.10;
	else if (shift=="lptclosure_xtrg_qcd_down") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_xtrg_taupt30to40)*0.90;
	else ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_xtrg_taupt30to40);
	if (shift=="lptclosure_xtrg_w_up") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_xtrg_taupt30to40)*1.10;
	else if (shift=="lptclosure_xtrg_w_down") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_xtrg_taupt30to40)*0.90;
	else ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_xtrg_taupt30to40);
	if (shift=="lptclosure_xtrg_tt_up") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_xtrg_taupt30to40)*1.10;
	else if (shift=="lptclosure_xtrg_tt_down") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_xtrg_taupt30to40)*0.90;
	else ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_xtrg_taupt30to40);
      }
      else if (pt<50){
	if (shift=="lptclosure_xtrg_qcd_up") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_xtrg_taupt40to50)*1.10;
	else if (shift=="lptclosure_xtrg_qcd_down") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_xtrg_taupt40to50)*0.90;
	else ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_xtrg_taupt40to50);
	if (shift=="lptclosure_xtrg_w_up") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_xtrg_taupt40to50)*1.10;
	else if (shift=="lptclosure_xtrg_w_down") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_xtrg_taupt40to50)*0.90;
	else ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_xtrg_taupt40to50);
	if (shift=="lptclosure_xtrg_tt_up") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_xtrg_taupt40to50)*1.10;
	else if (shift=="lptclosure_xtrg_tt_down") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_xtrg_taupt40to50)*0.90;
	else ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_xtrg_taupt40to50);
      }
      else{
	if (shift=="lptclosure_xtrg_qcd_up") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_xtrg_tauptgt50)*1.10;
	else if (shift=="lptclosure_xtrg_qcd_down") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_xtrg_tauptgt50)*0.90;
	else ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_xtrg_tauptgt50);
	if (shift=="lptclosure_xtrg_w_up") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_xtrg_tauptgt50)*1.10;
	else if (shift=="lptclosure_xtrg_w_down") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_xtrg_tauptgt50)*0.90;
	else ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_xtrg_tauptgt50);
	if (shift=="lptclosure_xtrg_tt_up") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_xtrg_tauptgt50)*1.10;
	else if (shift=="lptclosure_xtrg_tt_down") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_xtrg_tauptgt50)*0.90;
	else ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_xtrg_tauptgt50);
      }
    }
    else{
      if (pt<40){
	if (shift=="lptclosure_qcd_up") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_taupt30to40)*1.10;
	else if (shift=="lptclosure_qcd_down") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_taupt30to40)*0.90;
	else ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_taupt30to40);
	if (shift=="lptclosure_w_up") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_taupt30to40)*1.10;
	else if (shift=="lptclosure_w_down") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_taupt30to40)*0.90;
	else ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_taupt30to40);
	if (shift=="lptclosure_tt_up") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_taupt30to40)*1.10;
	else if (shift=="lptclosure_tt_down") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_taupt30to40)*0.90;
	else ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_taupt30to40);
      }
      else if (pt<50){
	if (shift=="lptclosure_qcd_up") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_taupt40to50)*1.10;
	else if (shift=="lptclosure_qcd_down") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_taupt40to50)*0.90;
	else ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_taupt40to50);
	if (shift=="lptclosure_w_up") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_taupt40to50)*1.10;
	else if (shift=="lptclosure_w_down") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_taupt40to50)*0.90;
	else ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_taupt40to50);
	if (shift=="lptclosure_tt_up") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_taupt40to50)*1.10;
	else if (shift=="lptclosure_tt_down") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_taupt40to50)*0.90;
	else ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_taupt40to50);
      }
      else{
	if (shift=="lptclosure_qcd_up") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_tauptgt50)*1.10;
	else if (shift=="lptclosure_qcd_down") ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_tauptgt50)*0.90;
	else ff_qcd = ff_qcd*get_lpt_closure(lpt, lptclosure_qcd_tauptgt50);
	if (shift=="lptclosure_w_up") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_tauptgt50)*1.10;
	else if (shift=="lptclosure_w_down") ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_tauptgt50)*0.90;
	else ff_w = ff_w*get_lpt_closure(lpt, lptclosure_w_tauptgt50);
	if (shift=="lptclosure_tt_up") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_tauptgt50)*1.10;
	else if (shift=="lptclosure_tt_down") ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_tauptgt50)*0.90;
	else ff_tt = ff_tt*get_lpt_closure(lpt, lptclosure_tt_tauptgt50);
      }
    }

    // MT and OSSS corrections
    if (shift=="mtclosure_w_up1") ff_w = ff_w*get_mt_closure(mt,mtclosure_w_up1);
    else if (shift=="mtclosure_w_up2") ff_w = ff_w*get_mt_closure(mt,mtclosure_w_up2);
    else if (shift=="mtclosure_w_down1") ff_w = ff_w*get_mt_closure(mt,mtclosure_w_down1);
    else if (shift=="mtclosure_w_down2") ff_w = ff_w*get_mt_closure(mt,mtclosure_w_down2);
    else ff_w = ff_w*get_mt_closure(mt,mtclosure_w);

    if (shift=="osssclosure_qcd_down") ff_qcd = ff_qcd*get_mvis_closure(3.0, osssclosure_qcd)*0.90;
    else if (shift=="osssclosure_qcd_up") ff_qcd = ff_qcd*get_mvis_closure(3.0, osssclosure_qcd)*1.10;
    else ff_qcd = ff_qcd*get_mvis_closure(3.0, osssclosure_qcd);
    /*if (shift=="osssclosure_qcd_up") ff_qcd = ff_qcd*1.3;
      else if (shift=="osssclosure_qcd_down") ff_qcd = ff_qcd;
      else ff_qcd = ff_qcd*1.15;*/ //FIXME

    float ff_cmb=frac_tt*ff_tt + frac_qcd*ff_qcd + frac_w*ff_w;
    return ff_cmb;
  } 

};
