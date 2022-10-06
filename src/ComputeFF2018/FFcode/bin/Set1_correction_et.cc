#include <TH2.h>
#include "ComputeFF2018/FFcode/interface/ApplyFF.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include "TMultiGraph.h"
#include <iostream>
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
#include "ComputeFF2018/FFcode/interface/ScaleFactor.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TKey.h"
#include "THashList.h"
#include "THStack.h"
#include "TPaveLabel.h"
#include "TFile.h"
#include "ComputeFF2018/FFcode/interface/myHelper.h"
#include "ComputeFF2018/FFcode/interface/et_Tree.h"
#include "ComputeFF2018/FFcode/interface/LumiReweightingStandAlone.h"
#include "TauAnalysisTools/TauTriggerSFs/interface/TauTriggerSFs2017.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"
#include "ComputeFF2018/FFcode/interface/SFtautrigger.h"
#include "TauPOG/TauIDSFs/interface/TauIDSFTool.h"

using namespace std;

int main(int argc, char** argv) {

    int do_control_plots=0;

    std::string input = *(argv + 1);
    std::string output = *(argv + 2);
    std::string sample = *(argv + 3);
    std::string name = *(argv + 4);
    std::string year = *(argv + 5);
    std::string apply_corrections = *(argv + 6);

    TFile *f_Double = new TFile(input.c_str());
    cout<<"XXXXXXXXXXXXX "<<input.c_str()<<" XXXXXXXXXXXX"<<endl;
    TTree *arbre = (TTree*) f_Double->Get("etau_tree");
    TH1F* nbevt = (TH1F*) f_Double->Get("nevents");
    float ngen = nbevt->GetBinContent(2);

    float xs=1.0; float weight=1.0; float luminosity=59740.0;
    if (year=="2017") luminosity=41529.0;
    if (year=="2016") luminosity=35922.0;
    if (sample=="DY" or sample=="ZL" or sample=="ZTT" or sample=="ZJ" or sample=="ZLL"){ xs=6225.42; weight=luminosity*xs/ngen;}
    if (sample=="DYlow"){ xs=21658.0; weight=luminosity*xs/ngen;}
    else if (sample=="TT" or sample=="TTT" or sample=="TTJ") {xs=831.76; weight=luminosity*xs/ngen;}
    else if (sample=="TTTo2L2Nu") {xs=88.29; weight=luminosity*xs/ngen;}
    else if (sample=="TTToSemiLeptonic") {xs=365.35; weight=luminosity*xs/ngen;}
    else if (sample=="TTToHadronic") {xs=377.96; weight=luminosity*xs/ngen;}
    else if (sample=="W") {xs=61526.7; weight=luminosity*xs/ngen;}
    else if (sample=="data_obs"){weight=1.0;}
    else if (sample=="embedded"){weight=1.0;}
    else if (sample=="WZ1L1Nu2Q") {xs=11.66; weight=luminosity*xs/ngen;}
    else if (sample=="WZ1L3Nu") {xs=3.293; weight=luminosity*xs/ngen;}
    else if (sample=="WZJets") {xs=5.26; weight=luminosity*xs/ngen;}
    else if (sample=="WZ3LNu") {xs=5.052; weight=luminosity*xs/ngen;}
    else if (sample=="WZ2L2Q") {xs=6.331; weight=luminosity*xs/ngen;}
    else if (sample=="WW1L1Nu2Q") {xs=49.997; weight=luminosity*xs/ngen;}
    else if (sample=="ZZ4L") {xs=1.325; weight=luminosity*xs/ngen;}
    else if (sample=="VV2L2Nu") {xs=13.97; weight=luminosity*xs/ngen;}
    else if (sample=="WW2L2Nu") {xs=11.08; weight=luminosity*xs/ngen;}
    else if (sample=="ZZ2L2Nu") {xs=0.6008; weight=luminosity*xs/ngen;}
    else if (sample=="ZZ2L2Q") {xs=3.688; weight=luminosity*xs/ngen;}
    else if (sample=="WW4Q") {xs=47.73; weight=luminosity*xs/ngen;}
    else if (sample=="WWLNuQQ") {xs=45.99; weight=luminosity*xs/ngen;}
    else if (sample=="ST_tW_antitop") {xs=35.85; weight=luminosity*xs/ngen;}
    else if (sample=="ST_tW_top") {xs=35.85; weight=luminosity*xs/ngen;}
    else if (sample=="ST_t_antitop") {xs=26.23; weight=luminosity*xs/ngen;}
    else if (sample=="ST_t_top") {xs=44.07; weight=luminosity*xs/ngen;}
    else if (sample=="ggH_htt125") {xs=48.58*0.0627; weight=luminosity*xs/ngen;}
    else if (sample=="qqH_htt125") {xs=3.782*0.0627; weight=luminosity*xs/ngen;}
    else if (sample=="WplusH125") {xs=0.840*0.0627; weight=luminosity*xs/ngen;}
    else if (sample=="WminusH125") {xs=0.5328*0.0627; weight=luminosity*xs/ngen;}
    else if (sample=="ZH125") {xs=0.8839*0.0627; weight=luminosity*xs/ngen;}
    else if (sample=="ZZ") {xs=16.523; weight=luminosity*xs/ngen;}
    else if (sample=="WZ") {xs=47.13; weight=luminosity*xs/ngen;}
    else if (sample=="WW") {xs=118.7; weight=luminosity*xs/ngen;}
    else if (sample=="WGLNu") {xs=489.0; weight=luminosity*xs/ngen;}
    else if (sample=="WGstarMuMu") {xs=2.793; weight=luminosity*xs/ngen;}
    else if (sample=="WGstarEE") {xs=3.526; weight=luminosity*xs/ngen;}
    else if (sample=="EWKWminus") {xs=23.24; weight=luminosity*xs/ngen;}
    else if (sample=="EWKWplus") {xs=29.59; weight=luminosity*xs/ngen;}
    else if (sample=="EWKZLL" or sample=="EWKZLL_TT" or sample=="EWKZLL_J" or sample=="EWKZLL_L" or sample=="EWKZLL_LL") {xs=4.321; weight=luminosity*xs/ngen;}
    else if (sample=="EWKZNuNu" or sample=="EWKZNuNu_TT" or sample=="EWKZNuNu_J" or sample=="EWKZNuNu_L" or sample=="EWKZNuNu_LL") {xs=10.66; weight=luminosity*xs/ngen;}

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(10);

    arbre->SetBranchAddress("Rivet_higgsPt", &Rivet_higgsPt);
    arbre->SetBranchAddress("Rivet_nJets30", &Rivet_nJets30);
    arbre->SetBranchAddress("Rivet_stage0_cat", &Rivet_stage0_cat);
    arbre->SetBranchAddress("Rivet_stage1p1_cat", &Rivet_stage1p1_cat);
    arbre->SetBranchAddress("Rivet_stage1_cat_pTjet30GeV", &Rivet_stage1_cat_pTjet30GeV);

    arbre->SetBranchAddress("passEle27", &passEle27);
    arbre->SetBranchAddress("passEle32", &passEle32);
    arbre->SetBranchAddress("passEle35", &passEle35);
    arbre->SetBranchAddress("passEle24Tau30", &passEle24Tau30);
    arbre->SetBranchAddress("passEle24HPSTau30", &passEle24HPSTau30);
    arbre->SetBranchAddress("passEle25", &passEle25);
    arbre->SetBranchAddress("matchEle25_1", &matchEle25_1);
    arbre->SetBranchAddress("filterEle25_1", &filterEle25_1);

    arbre->SetBranchAddress("matchEle27_1", &matchEle27_1);
    arbre->SetBranchAddress("matchEle32_1", &matchEle32_1);
    arbre->SetBranchAddress("matchEle35_1", &matchEle35_1);
    arbre->SetBranchAddress("matchEle24Tau30_1", &matchEle24Tau30_1);
    arbre->SetBranchAddress("matchEle24Tau30_2", &matchEle24Tau30_2);
    arbre->SetBranchAddress("matchEle24HPSTau30_1", &matchEle24HPSTau30_1);
    arbre->SetBranchAddress("matchEle24HPSTau30_2", &matchEle24HPSTau30_2);
    arbre->SetBranchAddress("filterEle27_1", &filterEle27_1);
    arbre->SetBranchAddress("filterEle32_1", &filterEle32_1);
    arbre->SetBranchAddress("filterEle35_1", &filterEle35_1);
    arbre->SetBranchAddress("filterEle24Tau30_1", &filterEle24Tau30_1);
    arbre->SetBranchAddress("filterEle24Tau30_2", &filterEle24Tau30_2);
    arbre->SetBranchAddress("filterEle24HPSTau30_1", &filterEle24HPSTau30_1);
    arbre->SetBranchAddress("filterEle24HPSTau30_2", &filterEle24HPSTau30_2);

    arbre->SetBranchAddress("run", &run);
    arbre->SetBranchAddress("lumi", &lumi);
    arbre->SetBranchAddress("evt", &evt);
    arbre->SetBranchAddress("npv", &npv);
    arbre->SetBranchAddress("pt_1", &pt_1);
    arbre->SetBranchAddress("phi_1", &phi_1);
    arbre->SetBranchAddress("eta_1", &eta_1);
    arbre->SetBranchAddress("iso_1", &iso_1);
    arbre->SetBranchAddress("isoDB_1", &isoDB_1);
    arbre->SetBranchAddress("m_1", &m_1);
    arbre->SetBranchAddress("q_1", &q_1);
    arbre->SetBranchAddress("nbtag", &nbtag);
    arbre->SetBranchAddress("nbtagL", &nbtagL);
    arbre->SetBranchAddress("bweight", &bweight);
    arbre->SetBranchAddress("prefiring_weight", &prefiring_weight);
    arbre->SetBranchAddress("q_2", &q_2);
    arbre->SetBranchAddress("pt_2", &pt_2);
    arbre->SetBranchAddress("eta_2", &eta_2);
    arbre->SetBranchAddress("m_2", &m_2);
    arbre->SetBranchAddress("phi_2", &phi_2);
    arbre->SetBranchAddress("met", &met);
    arbre->SetBranchAddress("m_sv", &m_sv);
    arbre->SetBranchAddress("m_sv_DOWN", &m_sv_DOWN);
    arbre->SetBranchAddress("m_sv_UP", &m_sv_UP);
    arbre->SetBranchAddress("m_sv_UESUp", &m_sv_UESUp);
    arbre->SetBranchAddress("m_sv_UESDown", &m_sv_UESDown);
    arbre->SetBranchAddress("met", &met);
    arbre->SetBranchAddress("metphi", &metphi);
    arbre->SetBranchAddress("met_UESUp", &met_UESUp);
    arbre->SetBranchAddress("metphi_UESDown", &metphi_UESDown);
    arbre->SetBranchAddress("met_UESDown", &met_UESDown);
    arbre->SetBranchAddress("metphi_UESUp", &metphi_UESUp);
    arbre->SetBranchAddress("met_responseUp", &met_responseUp);
    arbre->SetBranchAddress("metphi_responseUp", &metphi_responseUp);
    arbre->SetBranchAddress("met_responseDown", &met_responseDown);
    arbre->SetBranchAddress("metphi_responseDown", &metphi_responseDown);
    arbre->SetBranchAddress("met_resolutionUp", &met_resolutionUp);
    arbre->SetBranchAddress("metphi_resolutionUp", &metphi_resolutionUp);
    arbre->SetBranchAddress("met_resolutionDown", &met_resolutionDown);
    arbre->SetBranchAddress("metphi_resolutionDown", &metphi_resolutionDown);
    arbre->SetBranchAddress("metphi_JetEta0to3Down", &metphi_JetEta0to3Down);
    arbre->SetBranchAddress("met_JetEta0to3Down", &met_JetEta0to3Down);
    arbre->SetBranchAddress("metphi_JetEta0to3Up", &metphi_JetEta0to3Up);
    arbre->SetBranchAddress("met_JetEta0to3Up", &met_JetEta0to3Up);
    arbre->SetBranchAddress("metphi_JetEta0to5Down", &metphi_JetEta0to5Down);
    arbre->SetBranchAddress("met_JetEta0to5Down", &met_JetEta0to5Down);
    arbre->SetBranchAddress("metphi_JetEta0to5Up", &metphi_JetEta0to5Up);
    arbre->SetBranchAddress("met_JetEta0to5Up", &met_JetEta0to5Up);
    arbre->SetBranchAddress("metphi_JetEta3to5Down", &metphi_JetEta3to5Down);
    arbre->SetBranchAddress("met_JetEta3to5Down", &met_JetEta3to5Down);
    arbre->SetBranchAddress("metphi_JetEta3to5Up", &metphi_JetEta3to5Up);
    arbre->SetBranchAddress("met_JetEta3to5Up", &met_JetEta3to5Up);
    arbre->SetBranchAddress("metphi_JetRelativeSampleDown", &metphi_JetRelativeSampleDown);
    arbre->SetBranchAddress("met_JetRelativeSampleDown", &met_JetRelativeSampleDown);
    arbre->SetBranchAddress("metphi_JetRelativeSampleUp", &metphi_JetRelativeSampleUp);
    arbre->SetBranchAddress("met_JetRelativeSampleUp", &met_JetRelativeSampleUp);
    arbre->SetBranchAddress("metphi_JetRelativeBalDown", &metphi_JetRelativeBalDown);
    arbre->SetBranchAddress("met_JetRelativeBalDown", &met_JetRelativeBalDown);
    arbre->SetBranchAddress("metphi_JetRelativeBalUp", &metphi_JetRelativeBalUp);
    arbre->SetBranchAddress("met_JetRelativeBalUp", &met_JetRelativeBalUp);
    arbre->SetBranchAddress("metphi_JetEC2Down", &metphi_JetEC2Down);
    arbre->SetBranchAddress("met_JetEC2Down", &met_JetEC2Down);
    arbre->SetBranchAddress("metphi_JetEC2Up", &metphi_JetEC2Up);
    arbre->SetBranchAddress("met_JetEC2Up", &met_JetEC2Up);
    arbre->SetBranchAddress("njets_JetEC2Down", &njets_JetEC2Down);
    arbre->SetBranchAddress("njets_JetEC2Up", &njets_JetEC2Up);
    arbre->SetBranchAddress("mjj_JetEC2Down", &mjj_JetEC2Down);
    arbre->SetBranchAddress("mjj_JetEC2Up", &mjj_JetEC2Up);
    arbre->SetBranchAddress("pt_1_SigmaUp", &pt_1_SigmaUp);
    arbre->SetBranchAddress("pt_1_SigmaDown", &pt_1_SigmaDown);
    arbre->SetBranchAddress("pt_1_ScaleUp", &pt_1_ScaleUp);
    arbre->SetBranchAddress("pt_1_ScaleDown", &pt_1_ScaleDown);
    arbre->SetBranchAddress("m_sv_ResponseDown", &m_sv_ResponseDown);
    arbre->SetBranchAddress("m_sv_ResponseUp", &m_sv_ResponseUp);
    arbre->SetBranchAddress("m_sv_ResolutionDown", &m_sv_ResolutionDown);
    arbre->SetBranchAddress("m_sv_ResolutionUp", &m_sv_ResolutionUp);
    arbre->SetBranchAddress("m_sv_JetRelativeSampleUp", &m_sv_JetRelativeSampleUp);
    arbre->SetBranchAddress("m_sv_JetRelativeSampleDown", &m_sv_JetRelativeSampleDown);
    arbre->SetBranchAddress("m_sv_JetRelativeBalUp", &m_sv_JetRelativeBalUp);
    arbre->SetBranchAddress("m_sv_JetRelativeBalDown", &m_sv_JetRelativeBalDown);
    arbre->SetBranchAddress("m_sv_JetEta3to5Up", &m_sv_JetEta3to5Up);
    arbre->SetBranchAddress("m_sv_JetEta3to5Down", &m_sv_JetEta3to5Down);
    arbre->SetBranchAddress("m_sv_JetEta0to3Up", &m_sv_JetEta0to3Up);
    arbre->SetBranchAddress("m_sv_JetEta0to3Down", &m_sv_JetEta0to3Down);
    arbre->SetBranchAddress("m_sv_JetEta0to5Up", &m_sv_JetEta0to5Up);
    arbre->SetBranchAddress("m_sv_JetEta0to5Down", &m_sv_JetEta0to5Down);
    arbre->SetBranchAddress("m_sv_JetEC2Up", &m_sv_JetEC2Up);
    arbre->SetBranchAddress("m_sv_JetEC2Down", &m_sv_JetEC2Down);
    arbre->SetBranchAddress("m_sv_ESMEARUP", &m_sv_ESMEARUP);
    arbre->SetBranchAddress("m_sv_ESMEARDOWN", &m_sv_ESMEARDOWN);
    arbre->SetBranchAddress("m_sv_ESCALEUP", &m_sv_ESCALEUP);
    arbre->SetBranchAddress("m_sv_ESCALEDOWN", &m_sv_ESCALEDOWN);
    arbre->SetBranchAddress("njets", &njets);
    arbre->SetBranchAddress("njets_JetEta0to3Down", &njets_JetEta0to3Down);
    arbre->SetBranchAddress("njets_JetEta0to3Up", &njets_JetEta0to3Up);
    arbre->SetBranchAddress("njets_JetRelativeSampleDown", &njets_JetRelativeSampleDown);
    arbre->SetBranchAddress("njets_JetRelativeSampleUp", &njets_JetRelativeSampleUp);
    arbre->SetBranchAddress("njets_JetRelativeBalDown", &njets_JetRelativeBalDown);
    arbre->SetBranchAddress("njets_JetRelativeBalUp", &njets_JetRelativeBalUp);
    arbre->SetBranchAddress("njets_JetEta0to5Down", &njets_JetEta0to5Down);
    arbre->SetBranchAddress("njets_JetEta0to5Up", &njets_JetEta0to5Up);
    arbre->SetBranchAddress("njets_JetEta3to5Down", &njets_JetEta3to5Down);
    arbre->SetBranchAddress("njets_JetEta3to5Up", &njets_JetEta3to5Up);
    arbre->SetBranchAddress("jpt_1", &jpt_1);
    arbre->SetBranchAddress("jeta_1", &jeta_1);
    arbre->SetBranchAddress("jphi_1", &jphi_1);
    arbre->SetBranchAddress("jpt_2", &jpt_2);
    arbre->SetBranchAddress("jeta_2", &jeta_2);
    arbre->SetBranchAddress("jphi_2", &jphi_2);
    arbre->SetBranchAddress("bpt_1", &bpt_1);
    arbre->SetBranchAddress("beta_1", &beta_1);
    arbre->SetBranchAddress("bphi_1", &bphi_1);
    arbre->SetBranchAddress("bpt_2", &bpt_2);
    arbre->SetBranchAddress("beta_2", &beta_2);
    arbre->SetBranchAddress("bphi_2", &bphi_2);
    arbre->SetBranchAddress("bflavor_1", &bflavor_1);
    arbre->SetBranchAddress("bflavor_2", &bflavor_2);
    arbre->SetBranchAddress("mjj", &mjj);
    arbre->SetBranchAddress("mjj_JetEta0to3Down", &mjj_JetEta0to3Down);
    arbre->SetBranchAddress("mjj_JetEta0to3Up", &mjj_JetEta0to3Up);
    arbre->SetBranchAddress("mjj_JetEta3to5Down", &mjj_JetEta3to5Down);
    arbre->SetBranchAddress("mjj_JetEta3to5Up", &mjj_JetEta3to5Up);
    arbre->SetBranchAddress("mjj_JetEta0to5Down", &mjj_JetEta0to5Down);
    arbre->SetBranchAddress("mjj_JetEta0to5Up", &mjj_JetEta0to5Up);
    arbre->SetBranchAddress("mjj_JetRelativeBalDown", &mjj_JetRelativeBalDown);
    arbre->SetBranchAddress("mjj_JetRelativeBalUp", &mjj_JetRelativeBalUp);
    arbre->SetBranchAddress("mjj_JetRelativeSampleDown", &mjj_JetRelativeSampleDown);
    arbre->SetBranchAddress("mjj_JetRelativeSampleUp", &mjj_JetRelativeSampleUp);
    arbre->SetBranchAddress("genweight", &genweight);
    arbre->SetBranchAddress("byVLooseIsolationMVArun2v2DBoldDMwLT_2",&byVLooseIsolationMVArun2v2DBoldDMwLT_2);
    arbre->SetBranchAddress("byLooseIsolationMVArun2v2DBoldDMwLT_2",&byLooseIsolationMVArun2v2DBoldDMwLT_2);
    arbre->SetBranchAddress("byMediumIsolationMVArun2v2DBoldDMwLT_2",&byMediumIsolationMVArun2v2DBoldDMwLT_2);
    arbre->SetBranchAddress("byTightIsolationMVArun2v2DBoldDMwLT_2",&byTightIsolationMVArun2v2DBoldDMwLT_2);
    arbre->SetBranchAddress("byVTightIsolationMVArun2v2DBoldDMwLT_2",&byVTightIsolationMVArun2v2DBoldDMwLT_2);
    arbre->SetBranchAddress("byIsolationMVA3oldDMwLTraw_2",&byIsolationMVA3oldDMwLTraw_2);
    arbre->SetBranchAddress("l2_decayMode",&l2_decayMode);
    arbre->SetBranchAddress("againstElectronTightMVA6_2",&againstElectronTightMVA6_2);
    arbre->SetBranchAddress("againstElectronVTightMVA6_2",&againstElectronVTightMVA6_2);
    arbre->SetBranchAddress("againstElectronTightMVA62018_2",&againstElectronTightMVA62018_2);
    arbre->SetBranchAddress("againstElectronVTightMVA62018_2",&againstElectronVTightMVA62018_2);
    arbre->SetBranchAddress("againstMuonLoose3_2",&againstMuonLoose3_2);

    arbre->SetBranchAddress("byVVVLooseDeepVSjet_2",&byVVVLooseDeepVSjet_2);
    arbre->SetBranchAddress("byVVLooseDeepVSjet_2",&byVVLooseDeepVSjet_2);
    arbre->SetBranchAddress("byVLooseDeepVSjet_2",&byVLooseDeepVSjet_2);
    arbre->SetBranchAddress("byLooseDeepVSjet_2",&byLooseDeepVSjet_2);
    arbre->SetBranchAddress("byMediumDeepVSjet_2",&byMediumDeepVSjet_2);
    arbre->SetBranchAddress("byTightDeepVSjet_2",&byTightDeepVSjet_2);
    arbre->SetBranchAddress("byVTightDeepVSjet_2",&byVTightDeepVSjet_2);
    arbre->SetBranchAddress("byVVVLooseDeepVSmu_2",&byVVVLooseDeepVSmu_2);
    arbre->SetBranchAddress("byVVLooseDeepVSmu_2",&byVVLooseDeepVSmu_2);
    arbre->SetBranchAddress("byVLooseDeepVSmu_2",&byVLooseDeepVSmu_2);
    arbre->SetBranchAddress("byLooseDeepVSmu_2",&byLooseDeepVSmu_2);
    arbre->SetBranchAddress("byMediumDeepVSmu_2",&byMediumDeepVSmu_2);
    arbre->SetBranchAddress("byTightDeepVSmu_2",&byTightDeepVSmu_2);
    arbre->SetBranchAddress("byVTightDeepVSmu_2",&byVTightDeepVSmu_2);
    arbre->SetBranchAddress("byVVVLooseDeepVSe_2",&byVVVLooseDeepVSe_2);
    arbre->SetBranchAddress("byVVLooseDeepVSe_2",&byVVLooseDeepVSe_2);
    arbre->SetBranchAddress("byVLooseDeepVSe_2",&byVLooseDeepVSe_2);
    arbre->SetBranchAddress("byLooseDeepVSe_2",&byLooseDeepVSe_2);
    arbre->SetBranchAddress("byMediumDeepVSe_2",&byMediumDeepVSe_2);
    arbre->SetBranchAddress("byTightDeepVSe_2",&byTightDeepVSe_2);
    arbre->SetBranchAddress("byVTightDeepVSe_2",&byVTightDeepVSe_2);

    arbre->SetBranchAddress("gen_match_1",&gen_match_1);
    arbre->SetBranchAddress("gen_match_2",&gen_match_2);
    arbre->SetBranchAddress("npu",&npu);
    arbre->SetBranchAddress("genpT",&genpT);
    arbre->SetBranchAddress("genM",&genM);
    arbre->SetBranchAddress("pt_top1",&pt_top1);
    arbre->SetBranchAddress("pt_top2",&pt_top2);
    arbre->SetBranchAddress("numGenJets",&numGenJets);
    arbre->SetBranchAddress("Flag_BadChargedCandidateFilter",&Flag_BadChargedCandidateFilter);
    arbre->SetBranchAddress("Flag_BadPFMuonFilter",&Flag_BadPFMuonFilter);
    arbre->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter",&Flag_EcalDeadCellTriggerPrimitiveFilter);
    arbre->SetBranchAddress("Flag_HBHENoiseFilter",&Flag_HBHENoiseFilter);
    arbre->SetBranchAddress("Flag_HBHENoiseIsoFilter",&Flag_HBHENoiseIsoFilter);
    arbre->SetBranchAddress("Flag_eeBadScFilter",&Flag_eeBadScFilter);
    arbre->SetBranchAddress("Flag_goodVertices",&Flag_goodVertices);
    arbre->SetBranchAddress("Flag_globalSuperTightHalo2016Filter",&Flag_globalSuperTightHalo2016Filter);
    arbre->SetBranchAddress("Flag_ecalBadCalibFilter",&Flag_ecalBadCalibFilter);
    arbre->SetBranchAddress("Flag_duplicateMuons",&Flag_duplicateMuons);
    arbre->SetBranchAddress("Flag_badMuons",&Flag_badMuons);
    arbre->SetBranchAddress("Flag_ecalBadCalibReducedMINIAODFilter",&Flag_ecalBadCalibReducedMINIAODFilter);

    arbre->SetBranchAddress("matchEmbFilter_Ele24Tau30_1",&matchEmbFilter_Ele24Tau30_1);
    arbre->SetBranchAddress("matchEmbFilter_Ele27_1",&matchEmbFilter_Ele27_1);
    arbre->SetBranchAddress("matchEmbFilter_Ele32DoubleL1v1_1",&matchEmbFilter_Ele32DoubleL1v1_1);
    arbre->SetBranchAddress("matchEmbFilter_Ele32DoubleL1v2_1",&matchEmbFilter_Ele32DoubleL1v2_1);
    arbre->SetBranchAddress("matchEmbFilter_Ele32_1",&matchEmbFilter_Ele32_1);
    arbre->SetBranchAddress("matchEmbFilter_Ele35_1",&matchEmbFilter_Ele35_1);
    arbre->SetBranchAddress("matchEmbFilter_Ele24Tau30_2",&matchEmbFilter_Ele24Tau30_2);

   int nbhist=1;

   float bins_mtt0[] = {0,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350,400};
   int  binnum_mtt0 = sizeof(bins_mtt0)/sizeof(Float_t) - 1;
   float bins_taupt0[] = {30,35,40,45,50,55,60,65,70,75,80,120};
   int  binnum_taupt0 = sizeof(bins_taupt0)/sizeof(Float_t) - 1;
   float bins_pth0[] = {0,25,50,75,100,125,150,175,200,300,500};
   int  binnum_pth0 = sizeof(bins_pth0)/sizeof(Float_t) - 1;
   float bins_mjj0[] = {0,100,200,300,400,500,600,700,1000,1500};
   int  binnum_mjj0 = sizeof(bins_mjj0)/sizeof(Float_t) - 1;
   float bins_mt[] = {0,10,20,30,40,50,60,70,80,90,100,110,120};
   int  binnum_mt = sizeof(bins_mt)/sizeof(Float_t) - 1;

   float bins_cat0jetlow0[] = {30,40,50,9000};
   int  binnum_cat0jetlow0 = sizeof(bins_cat0jetlow0)/sizeof(Float_t) - 1;
   float bins_cat0jetlow1[] = {50,70,90,110,130,150,9000};
   int  binnum_cat0jetlow1 = sizeof(bins_cat0jetlow1)/sizeof(Float_t) - 1;

   float bins_cat0jethigh0[] = {30,40,50,60,9000};
   int  binnum_cat0jethigh0 = sizeof(bins_cat0jethigh0)/sizeof(Float_t) - 1;
   float bins_cat0jethigh1[] = {50,70,90,110,130,150,9000};
   int  binnum_cat0jethigh1 = sizeof(bins_cat0jethigh1)/sizeof(Float_t) - 1;

   float bins_catboosted10[] = {0,60,120,9000};
   int  binnum_catboosted10 = sizeof(bins_catboosted10)/sizeof(Float_t) - 1;
   float bins_catboosted11[] = {50,70,90,110,130,150,170,210,250,9000};
   int  binnum_catboosted11 = sizeof(bins_catboosted11)/sizeof(Float_t) - 1;

   float bins_catboosted20[] = {0,60,120,200,9000};
   int  binnum_catboosted20 = sizeof(bins_catboosted20)/sizeof(Float_t) - 1;
   float bins_catboosted21[] = {50,70,90,110,130,150,170,210,250,9000};
   int  binnum_catboosted21 = sizeof(bins_catboosted21)/sizeof(Float_t) - 1;

   float bins_catvbflow0[] = {350,700,1000,1500,9000};
   int  binnum_catvbflow0 = sizeof(bins_catvbflow0)/sizeof(Float_t) - 1;
   float bins_catvbflow1[] = {50,90,130,170,250,9000};
   int  binnum_catvbflow1 = sizeof(bins_catvbflow1)/sizeof(Float_t) - 1;

   float bins_catvbfhigh0[] = {350,9000};
   int  binnum_catvbfhigh0 = sizeof(bins_catvbfhigh0)/sizeof(Float_t) - 1;
   float bins_catvbfhigh1[] = {50,100,150,200,250,9000};
   int  binnum_catvbfhigh1 = sizeof(bins_catvbfhigh1)/sizeof(Float_t) - 1;

   TH1F* h0LT_qcd_mvis_iso = new TH1F ("h0LT_qcd_mvis_iso","h0LT_qcd_mvis_iso",binnum_mtt0,bins_mtt0); h0LT_qcd_mvis_iso->Sumw2();
   TH1F* h0LT_qcd_mvis_anti = new TH1F ("h0LT_qcd_mvis_anti","h0LT_qcd_mvis_anti",binnum_mtt0,bins_mtt0); h0LT_qcd_mvis_anti->Sumw2();
   TH1F* h0LT_w_mvis_iso = new TH1F ("h0LT_w_mvis_iso","h0LT_w_mvis_iso",binnum_mtt0,bins_mtt0); h0LT_w_mvis_iso->Sumw2();
   TH1F* h0LT_w_mvis_anti = new TH1F ("h0LT_w_mvis_anti","h0LT_w_mvis_anti",binnum_mtt0,bins_mtt0); h0LT_w_mvis_anti->Sumw2();
   TH1F* h0J_qcd_mvis_iso = new TH1F ("h0J_qcd_mvis_iso","h0J_qcd_mvis_iso",binnum_mtt0,bins_mtt0); h0J_qcd_mvis_iso->Sumw2();
   TH1F* h0J_qcd_mvis_anti = new TH1F ("h0J_qcd_mvis_anti","h0J_qcd_mvis_anti",binnum_mtt0,bins_mtt0); h0J_qcd_mvis_anti->Sumw2();
   TH1F* h0J_w_mvis_iso = new TH1F ("h0J_w_mvis_iso","h0J_w_mvis_iso",binnum_mtt0,bins_mtt0); h0J_w_mvis_iso->Sumw2();
   TH1F* h0J_w_mvis_anti = new TH1F ("h0J_w_mvis_anti","h0J_w_mvis_anti",binnum_mtt0,bins_mtt0); h0J_w_mvis_anti->Sumw2();

   TH1F* h1LT_qcd_mvis_iso = new TH1F ("h1LT_qcd_mvis_iso","h1LT_qcd_mvis_iso",binnum_mtt0,bins_mtt0); h1LT_qcd_mvis_iso->Sumw2();
   TH1F* h1LT_qcd_mvis_anti = new TH1F ("h1LT_qcd_mvis_anti","h1LT_qcd_mvis_anti",binnum_mtt0,bins_mtt0); h1LT_qcd_mvis_anti->Sumw2();
   TH1F* h1LT_w_mvis_iso = new TH1F ("h1LT_w_mvis_iso","h1LT_w_mvis_iso",binnum_mtt0,bins_mtt0); h1LT_w_mvis_iso->Sumw2();
   TH1F* h1LT_w_mvis_anti = new TH1F ("h1LT_w_mvis_anti","h1LT_w_mvis_anti",binnum_mtt0,bins_mtt0); h1LT_w_mvis_anti->Sumw2();
   TH1F* h1J_qcd_mvis_iso = new TH1F ("h1J_qcd_mvis_iso","h1J_qcd_mvis_iso",binnum_mtt0,bins_mtt0); h1J_qcd_mvis_iso->Sumw2();
   TH1F* h1J_qcd_mvis_anti = new TH1F ("h1J_qcd_mvis_anti","h1J_qcd_mvis_anti",binnum_mtt0,bins_mtt0); h1J_qcd_mvis_anti->Sumw2();
   TH1F* h1J_w_mvis_iso = new TH1F ("h1J_w_mvis_iso","h1J_w_mvis_iso",binnum_mtt0,bins_mtt0); h1J_w_mvis_iso->Sumw2();
   TH1F* h1J_w_mvis_anti = new TH1F ("h1J_w_mvis_anti","h1J_w_mvis_anti",binnum_mtt0,bins_mtt0); h1J_w_mvis_anti->Sumw2();

   TH1F* h2LT_qcd_mvis_iso = new TH1F ("h2LT_qcd_mvis_iso","h2LT_qcd_mvis_iso",binnum_mtt0,bins_mtt0); h2LT_qcd_mvis_iso->Sumw2();
   TH1F* h2LT_qcd_mvis_anti = new TH1F ("h2LT_qcd_mvis_anti","h2LT_qcd_mvis_anti",binnum_mtt0,bins_mtt0); h2LT_qcd_mvis_anti->Sumw2();
   TH1F* h2LT_w_mvis_iso = new TH1F ("h2LT_w_mvis_iso","h2LT_w_mvis_iso",binnum_mtt0,bins_mtt0); h2LT_w_mvis_iso->Sumw2();
   TH1F* h2LT_w_mvis_anti = new TH1F ("h2LT_w_mvis_anti","h2LT_w_mvis_anti",binnum_mtt0,bins_mtt0); h2LT_w_mvis_anti->Sumw2();
   TH1F* h2J_qcd_mvis_iso = new TH1F ("h2J_qcd_mvis_iso","h2J_qcd_mvis_iso",binnum_mtt0,bins_mtt0); h2J_qcd_mvis_iso->Sumw2();
   TH1F* h2J_qcd_mvis_anti = new TH1F ("h2J_qcd_mvis_anti","h2J_qcd_mvis_anti",binnum_mtt0,bins_mtt0); h2J_qcd_mvis_anti->Sumw2();
   TH1F* h2J_w_mvis_iso = new TH1F ("h2J_w_mvis_iso","h2J_w_mvis_iso",binnum_mtt0,bins_mtt0); h2J_w_mvis_iso->Sumw2();
   TH1F* h2J_w_mvis_anti = new TH1F ("h2J_w_mvis_anti","h2J_w_mvis_anti",binnum_mtt0,bins_mtt0); h2J_w_mvis_anti->Sumw2();


   TH1F* h0LT_qcd_iso = new TH1F ("h0LT_qcd_iso","h0LT_qcd_iso",binnum_mtt0,bins_mtt0); h0LT_qcd_iso->Sumw2();
   TH1F* h0LT_qcd_anti = new TH1F ("h0LT_qcd_anti","h0LT_qcd_anti",binnum_mtt0,bins_mtt0); h0LT_qcd_anti->Sumw2();
   TH1F* h0LT_qcd_taupt_iso = new TH1F ("h0LT_qcd_taupt_iso","h0LT_qcd_taupt_iso",binnum_taupt0,bins_taupt0); h0LT_qcd_taupt_iso->Sumw2();
   TH1F* h0LT_qcd_taupt_anti = new TH1F ("h0LT_qcd_taupt_anti","h0LT_qcd_taupt_anti",binnum_taupt0,bins_taupt0); h0LT_qcd_taupt_anti->Sumw2();
   TH1F* h0LT_qcd_pth_iso = new TH1F ("h0LT_qcd_pth_iso","h0LT_qcd_pth_iso",binnum_pth0,bins_pth0); h0LT_qcd_pth_iso->Sumw2();
   TH1F* h0LT_qcd_pth_anti = new TH1F ("h0LT_qcd_pth_anti","h0LT_qcd_pth_anti",binnum_pth0,bins_pth0); h0LT_qcd_pth_anti->Sumw2();
   TH1F* h0LT_qcd_mjj_iso = new TH1F ("h0LT_qcd_mjj_iso","h0LT_qcd_mjj_iso",binnum_mjj0,bins_mjj0); h0LT_qcd_mjj_iso->Sumw2();
   TH1F* h0LT_qcd_mjj_anti = new TH1F ("h0LT_qcd_mjj_anti","h0LT_qcd_mjj_anti",binnum_mjj0,bins_mjj0); h0LT_qcd_mjj_anti->Sumw2();
   TH1F* h0SSlooseLT_qcd_iso = new TH1F ("h0SSlooseLT_qcd_iso","h0SSlooseLT_qcd_iso",binnum_mtt0,bins_mtt0); h0SSlooseLT_qcd_iso->Sumw2();
   TH1F* h0SSlooseLT_qcd_anti = new TH1F ("h0SSlooseLT_qcd_anti","h0SSlooseLT_qcd_anti",binnum_mtt0,bins_mtt0); h0SSlooseLT_qcd_anti->Sumw2();
   TH1F* h0LT_w_iso = new TH1F ("h0LT_w_iso","h0LT_w_iso",binnum_mtt0,bins_mtt0); h0LT_w_iso->Sumw2();
   TH1F* h0LT_w_anti = new TH1F ("h0LT_w_anti","h0LT_w_anti",binnum_mtt0,bins_mtt0); h0LT_w_anti->Sumw2();
   TH1F* h0LT_w_taupt_iso = new TH1F ("h0LT_w_taupt_iso","h0LT_w_taupt_iso",binnum_taupt0,bins_taupt0); h0LT_w_taupt_iso->Sumw2();
   TH1F* h0LT_w_taupt_anti = new TH1F ("h0LT_w_taupt_anti","h0LT_w_taupt_anti",binnum_taupt0,bins_taupt0); h0LT_w_taupt_anti->Sumw2();
   TH1F* h0LT_w_pth_iso = new TH1F ("h0LT_w_pth_iso","h0LT_w_pth_iso",binnum_pth0,bins_pth0); h0LT_w_pth_iso->Sumw2();
   TH1F* h0LT_w_pth_anti = new TH1F ("h0LT_w_pth_anti","h0LT_w_pth_anti",binnum_pth0,bins_pth0); h0LT_w_pth_anti->Sumw2();
   TH1F* h0LT_w_mjj_iso = new TH1F ("h0LT_w_mjj_iso","h0LT_w_mjj_iso",binnum_mjj0,bins_mjj0); h0LT_w_mjj_iso->Sumw2();
   TH1F* h0LT_w_mjj_anti = new TH1F ("h0LT_w_mjj_anti","h0LT_w_mjj_anti",binnum_mjj0,bins_mjj0); h0LT_w_mjj_anti->Sumw2();
   TH1F* h0LT_tt_iso = new TH1F ("h0LT_tt_iso","h0LT_tt_iso",binnum_mtt0,bins_mtt0); h0LT_tt_iso->Sumw2();
   TH1F* h0LT_tt_anti = new TH1F ("h0LT_tt_anti","h0LT_tt_anti",binnum_mtt0,bins_mtt0); h0LT_tt_anti->Sumw2();

   TH1F* h0J_qcd_iso = new TH1F ("h0J_qcd_iso","h0J_qcd_iso",binnum_mtt0,bins_mtt0); h0J_qcd_iso->Sumw2();
   TH1F* h0J_qcd_anti = new TH1F ("h0J_qcd_anti","h0J_qcd_anti",binnum_mtt0,bins_mtt0); h0J_qcd_anti->Sumw2();
   TH1F* h0J_qcd_taupt_iso = new TH1F ("h0J_qcd_taupt_iso","h0J_qcd_taupt_iso",binnum_taupt0,bins_taupt0); h0J_qcd_taupt_iso->Sumw2();
   TH1F* h0J_qcd_taupt_anti = new TH1F ("h0J_qcd_taupt_anti","h0J_qcd_taupt_anti",binnum_taupt0,bins_taupt0); h0J_qcd_taupt_anti->Sumw2();
   TH1F* h0J_qcd_pth_iso = new TH1F ("h0J_qcd_pth_iso","h0J_qcd_pth_iso",binnum_pth0,bins_pth0); h0J_qcd_pth_iso->Sumw2();
   TH1F* h0J_qcd_pth_anti = new TH1F ("h0J_qcd_pth_anti","h0J_qcd_pth_anti",binnum_pth0,bins_pth0); h0J_qcd_pth_anti->Sumw2();
   TH1F* h0J_qcd_mjj_iso = new TH1F ("h0J_qcd_mjj_iso","h0J_qcd_mjj_iso",binnum_mjj0,bins_mjj0); h0J_qcd_mjj_iso->Sumw2();
   TH1F* h0J_qcd_mjj_anti = new TH1F ("h0J_qcd_mjj_anti","h0J_qcd_mjj_anti",binnum_mjj0,bins_mjj0); h0J_qcd_mjj_anti->Sumw2();
   TH1F* h0SSlooseJ_qcd_iso = new TH1F ("h0SSlooseJ_qcd_iso","h0SSlooseJ_qcd_iso",binnum_mtt0,bins_mtt0); h0SSlooseJ_qcd_iso->Sumw2();
   TH1F* h0SSlooseJ_qcd_anti = new TH1F ("h0SSlooseJ_qcd_anti","h0SSlooseJ_qcd_anti",binnum_mtt0,bins_mtt0); h0SSlooseJ_qcd_anti->Sumw2();
   TH1F* h0J_w_iso = new TH1F ("h0J_w_iso","h0J_w_iso",binnum_mtt0,bins_mtt0); h0J_w_iso->Sumw2();
   TH1F* h0J_w_anti = new TH1F ("h0J_w_anti","h0J_w_anti",binnum_mtt0,bins_mtt0); h0J_w_anti->Sumw2();
   TH1F* h0J_w_taupt_iso = new TH1F ("h0J_w_taupt_iso","h0J_w_taupt_iso",binnum_taupt0,bins_taupt0); h0J_w_taupt_iso->Sumw2();
   TH1F* h0J_w_taupt_anti = new TH1F ("h0J_w_taupt_anti","h0J_w_taupt_anti",binnum_taupt0,bins_taupt0); h0J_w_taupt_anti->Sumw2();
   TH1F* h0J_w_pth_iso = new TH1F ("h0J_w_pth_iso","h0J_w_pth_iso",binnum_pth0,bins_pth0); h0J_w_pth_iso->Sumw2();
   TH1F* h0J_w_pth_anti = new TH1F ("h0J_w_pth_anti","h0J_w_pth_anti",binnum_pth0,bins_pth0); h0J_w_pth_anti->Sumw2();
   TH1F* h0J_w_mjj_iso = new TH1F ("h0J_w_mjj_iso","h0J_w_mjj_iso",binnum_mjj0,bins_mjj0); h0J_w_mjj_iso->Sumw2();
   TH1F* h0J_w_mjj_anti = new TH1F ("h0J_w_mjj_anti","h0J_w_mjj_anti",binnum_mjj0,bins_mjj0); h0J_w_mjj_anti->Sumw2();
   TH1F* h0J_tt_iso = new TH1F ("h0J_tt_iso","h0J_tt_iso",binnum_mtt0,bins_mtt0); h0J_tt_iso->Sumw2();
   TH1F* h0J_tt_anti = new TH1F ("h0J_tt_anti","h0J_tt_anti",binnum_mtt0,bins_mtt0); h0J_tt_anti->Sumw2();

   TH1F* hmt_w_iso = new TH1F ("hmt_w_iso","hmt_w_iso",binnum_mt,bins_mt); hmt_w_iso->Sumw2();
   TH1F* hmt_w_anti = new TH1F ("hmt_w_anti","hmt_w_anti",binnum_mt,bins_mt); hmt_w_anti->Sumw2();

   TH2F* h0LT_qcd_cat0jetlow_iso = new TH2F ("h0LT_qcd_cat0jetlow_iso","h0LT_qcd_cat0jetlow_iso",binnum_cat0jetlow0,bins_cat0jetlow0,binnum_cat0jetlow1,bins_cat0jetlow1); h0LT_qcd_cat0jetlow_iso->Sumw2();
   TH2F* h0LT_qcd_cat0jetlow_anti = new TH2F ("h0LT_qcd_cat0jetlow_anti","h0LT_qcd_cat0jetlow_anti",binnum_cat0jetlow0,bins_cat0jetlow0,binnum_cat0jetlow1,bins_cat0jetlow1); h0LT_qcd_cat0jetlow_anti->Sumw2();
   TH2F* h0LT_w_cat0jetlow_iso = new TH2F ("h0LT_w_cat0jetlow_iso","h0LT_w_cat0jetlow_iso",binnum_cat0jetlow0,bins_cat0jetlow0,binnum_cat0jetlow1,bins_cat0jetlow1); h0LT_w_cat0jetlow_iso->Sumw2();
   TH2F* h0LT_w_cat0jetlow_anti = new TH2F ("h0LT_w_cat0jetlow_anti","h0LT_w_cat0jetlow_anti",binnum_cat0jetlow0,bins_cat0jetlow0,binnum_cat0jetlow1,bins_cat0jetlow1); h0LT_w_cat0jetlow_anti->Sumw2();
   TH2F* h0J_qcd_cat0jetlow_iso = new TH2F ("h0J_qcd_cat0jetlow_iso","h0J_qcd_cat0jetlow_iso",binnum_cat0jetlow0,bins_cat0jetlow0,binnum_cat0jetlow1,bins_cat0jetlow1); h0J_qcd_cat0jetlow_iso->Sumw2();
   TH2F* h0J_qcd_cat0jetlow_anti = new TH2F ("h0J_qcd_cat0jetlow_anti","h0J_qcd_cat0jetlow_anti",binnum_cat0jetlow0,bins_cat0jetlow0,binnum_cat0jetlow1,bins_cat0jetlow1); h0J_qcd_cat0jetlow_anti->Sumw2();
   TH2F* h0J_w_cat0jetlow_iso = new TH2F ("h0J_w_cat0jetlow_iso","h0J_w_cat0jetlow_iso",binnum_cat0jetlow0,bins_cat0jetlow0,binnum_cat0jetlow1,bins_cat0jetlow1); h0J_w_cat0jetlow_iso->Sumw2();
   TH2F* h0J_w_cat0jetlow_anti = new TH2F ("h0J_w_cat0jetlow_anti","h0J_w_cat0jetlow_anti",binnum_cat0jetlow0,bins_cat0jetlow0,binnum_cat0jetlow1,bins_cat0jetlow1); h0J_w_cat0jetlow_anti->Sumw2();

   TH2F* h0LT_qcd_cat0jethigh_iso = new TH2F ("h0LT_qcd_cat0jethigh_iso","h0LT_qcd_cat0jethigh_iso",binnum_cat0jethigh0,bins_cat0jethigh0,binnum_cat0jethigh1,bins_cat0jethigh1); h0LT_qcd_cat0jethigh_iso->Sumw2();
   TH2F* h0LT_qcd_cat0jethigh_anti = new TH2F ("h0LT_qcd_cat0jethigh_anti","h0LT_qcd_cat0jethigh_anti",binnum_cat0jethigh0,bins_cat0jethigh0,binnum_cat0jethigh1,bins_cat0jethigh1); h0LT_qcd_cat0jethigh_anti->Sumw2();
   TH2F* h0LT_w_cat0jethigh_iso = new TH2F ("h0LT_w_cat0jethigh_iso","h0LT_w_cat0jethigh_iso",binnum_cat0jethigh0,bins_cat0jethigh0,binnum_cat0jethigh1,bins_cat0jethigh1); h0LT_w_cat0jethigh_iso->Sumw2();
   TH2F* h0LT_w_cat0jethigh_anti = new TH2F ("h0LT_w_cat0jethigh_anti","h0LT_w_cat0jethigh_anti",binnum_cat0jethigh0,bins_cat0jethigh0,binnum_cat0jethigh1,bins_cat0jethigh1); h0LT_w_cat0jethigh_anti->Sumw2();
   TH2F* h0J_qcd_cat0jethigh_iso = new TH2F ("h0J_qcd_cat0jethigh_iso","h0J_qcd_cat0jethigh_iso",binnum_cat0jethigh0,bins_cat0jethigh0,binnum_cat0jethigh1,bins_cat0jethigh1); h0J_qcd_cat0jethigh_iso->Sumw2();
   TH2F* h0J_qcd_cat0jethigh_anti = new TH2F ("h0J_qcd_cat0jethigh_anti","h0J_qcd_cat0jethigh_anti",binnum_cat0jethigh0,bins_cat0jethigh0,binnum_cat0jethigh1,bins_cat0jethigh1); h0J_qcd_cat0jethigh_anti->Sumw2();
   TH2F* h0J_w_cat0jethigh_iso = new TH2F ("h0J_w_cat0jethigh_iso","h0J_w_cat0jethigh_iso",binnum_cat0jethigh0,bins_cat0jethigh0,binnum_cat0jethigh1,bins_cat0jethigh1); h0J_w_cat0jethigh_iso->Sumw2();
   TH2F* h0J_w_cat0jethigh_anti = new TH2F ("h0J_w_cat0jethigh_anti","h0J_w_cat0jethigh_anti",binnum_cat0jethigh0,bins_cat0jethigh0,binnum_cat0jethigh1,bins_cat0jethigh1); h0J_w_cat0jethigh_anti->Sumw2();

   TH2F* h0LT_qcd_catboosted1_iso = new TH2F ("h0LT_qcd_catboosted1_iso","h0LT_qcd_catboosted1_iso",binnum_catboosted10,bins_catboosted10,binnum_catboosted11,bins_catboosted11); h0LT_qcd_catboosted1_iso->Sumw2();
   TH2F* h0LT_qcd_catboosted1_anti = new TH2F ("h0LT_qcd_catboosted1_anti","h0LT_qcd_catboosted1_anti",binnum_catboosted10,bins_catboosted10,binnum_catboosted11,bins_catboosted11); h0LT_qcd_catboosted1_anti->Sumw2();
   TH2F* h0LT_w_catboosted1_iso = new TH2F ("h0LT_w_catboosted1_iso","h0LT_w_catboosted1_iso",binnum_catboosted10,bins_catboosted10,binnum_catboosted11,bins_catboosted11); h0LT_w_catboosted1_iso->Sumw2();
   TH2F* h0LT_w_catboosted1_anti = new TH2F ("h0LT_w_catboosted1_anti","h0LT_w_catboosted1_anti",binnum_catboosted10,bins_catboosted10,binnum_catboosted11,bins_catboosted11); h0LT_w_catboosted1_anti->Sumw2();
   TH2F* h0J_qcd_catboosted1_iso = new TH2F ("h0J_qcd_catboosted1_iso","h0J_qcd_catboosted1_iso",binnum_catboosted10,bins_catboosted10,binnum_catboosted11,bins_catboosted11); h0J_qcd_catboosted1_iso->Sumw2();
   TH2F* h0J_qcd_catboosted1_anti = new TH2F ("h0J_qcd_catboosted1_anti","h0J_qcd_catboosted1_anti",binnum_catboosted10,bins_catboosted10,binnum_catboosted11,bins_catboosted11); h0J_qcd_catboosted1_anti->Sumw2();
   TH2F* h0J_w_catboosted1_iso = new TH2F ("h0J_w_catboosted1_iso","h0J_w_catboosted1_iso",binnum_catboosted10,bins_catboosted10,binnum_catboosted11,bins_catboosted11); h0J_w_catboosted1_iso->Sumw2();
   TH2F* h0J_w_catboosted1_anti = new TH2F ("h0J_w_catboosted1_anti","h0J_w_catboosted1_anti",binnum_catboosted10,bins_catboosted10,binnum_catboosted11,bins_catboosted11); h0J_w_catboosted1_anti->Sumw2();

   TH2F* h0LT_qcd_catboosted2_iso = new TH2F ("h0LT_qcd_catboosted2_iso","h0LT_qcd_catboosted2_iso",binnum_catboosted20,bins_catboosted20,binnum_catboosted21,bins_catboosted21); h0LT_qcd_catboosted2_iso->Sumw2();
   TH2F* h0LT_qcd_catboosted2_anti = new TH2F ("h0LT_qcd_catboosted2_anti","h0LT_qcd_catboosted2_anti",binnum_catboosted20,bins_catboosted20,binnum_catboosted21,bins_catboosted21); h0LT_qcd_catboosted2_anti->Sumw2();
   TH2F* h0LT_w_catboosted2_iso = new TH2F ("h0LT_w_catboosted2_iso","h0LT_w_catboosted2_iso",binnum_catboosted20,bins_catboosted20,binnum_catboosted21,bins_catboosted21); h0LT_w_catboosted2_iso->Sumw2();
   TH2F* h0LT_w_catboosted2_anti = new TH2F ("h0LT_w_catboosted2_anti","h0LT_w_catboosted2_anti",binnum_catboosted20,bins_catboosted20,binnum_catboosted21,bins_catboosted21); h0LT_w_catboosted2_anti->Sumw2();
   TH2F* h0J_qcd_catboosted2_iso = new TH2F ("h0J_qcd_catboosted2_iso","h0J_qcd_catboosted2_iso",binnum_catboosted20,bins_catboosted20,binnum_catboosted21,bins_catboosted21); h0J_qcd_catboosted2_iso->Sumw2();
   TH2F* h0J_qcd_catboosted2_anti = new TH2F ("h0J_qcd_catboosted2_anti","h0J_qcd_catboosted2_anti",binnum_catboosted20,bins_catboosted20,binnum_catboosted21,bins_catboosted21); h0J_qcd_catboosted2_anti->Sumw2();
   TH2F* h0J_w_catboosted2_iso = new TH2F ("h0J_w_catboosted2_iso","h0J_w_catboosted2_iso",binnum_catboosted20,bins_catboosted20,binnum_catboosted21,bins_catboosted21); h0J_w_catboosted2_iso->Sumw2();
   TH2F* h0J_w_catboosted2_anti = new TH2F ("h0J_w_catboosted2_anti","h0J_w_catboosted2_anti",binnum_catboosted20,bins_catboosted20,binnum_catboosted21,bins_catboosted21); h0J_w_catboosted2_anti->Sumw2();

   TH2F* h0LT_qcd_catvbflow_iso = new TH2F ("h0LT_qcd_catvbflow_iso","h0LT_qcd_catvbflow_iso",binnum_catvbflow0,bins_catvbflow0,binnum_catvbflow1,bins_catvbflow1); h0LT_qcd_catvbflow_iso->Sumw2();
   TH2F* h0LT_qcd_catvbflow_anti = new TH2F ("h0LT_qcd_catvbflow_anti","h0LT_qcd_catvbflow_anti",binnum_catvbflow0,bins_catvbflow0,binnum_catvbflow1,bins_catvbflow1); h0LT_qcd_catvbflow_anti->Sumw2();
   TH2F* h0LT_w_catvbflow_iso = new TH2F ("h0LT_w_catvbflow_iso","h0LT_w_catvbflow_iso",binnum_catvbflow0,bins_catvbflow0,binnum_catvbflow1,bins_catvbflow1); h0LT_w_catvbflow_iso->Sumw2();
   TH2F* h0LT_w_catvbflow_anti = new TH2F ("h0LT_w_catvbflow_anti","h0LT_w_catvbflow_anti",binnum_catvbflow0,bins_catvbflow0,binnum_catvbflow1,bins_catvbflow1); h0LT_w_catvbflow_anti->Sumw2();
   TH2F* h0J_qcd_catvbflow_iso = new TH2F ("h0J_qcd_catvbflow_iso","h0J_qcd_catvbflow_iso",binnum_catvbflow0,bins_catvbflow0,binnum_catvbflow1,bins_catvbflow1); h0J_qcd_catvbflow_iso->Sumw2();
   TH2F* h0J_qcd_catvbflow_anti = new TH2F ("h0J_qcd_catvbflow_anti","h0J_qcd_catvbflow_anti",binnum_catvbflow0,bins_catvbflow0,binnum_catvbflow1,bins_catvbflow1); h0J_qcd_catvbflow_anti->Sumw2();
   TH2F* h0J_w_catvbflow_iso = new TH2F ("h0J_w_catvbflow_iso","h0J_w_catvbflow_iso",binnum_catvbflow0,bins_catvbflow0,binnum_catvbflow1,bins_catvbflow1); h0J_w_catvbflow_iso->Sumw2();
   TH2F* h0J_w_catvbflow_anti = new TH2F ("h0J_w_catvbflow_anti","h0J_w_catvbflow_anti",binnum_catvbflow0,bins_catvbflow0,binnum_catvbflow1,bins_catvbflow1); h0J_w_catvbflow_anti->Sumw2();

   TH2F* h0LT_qcd_catvbfhigh_iso = new TH2F ("h0LT_qcd_catvbfhigh_iso","h0LT_qcd_catvbfhigh_iso",binnum_catvbfhigh0,bins_catvbfhigh0,binnum_catvbfhigh1,bins_catvbfhigh1); h0LT_qcd_catvbfhigh_iso->Sumw2();
   TH2F* h0LT_qcd_catvbfhigh_anti = new TH2F ("h0LT_qcd_catvbfhigh_anti","h0LT_qcd_catvbfhigh_anti",binnum_catvbfhigh0,bins_catvbfhigh0,binnum_catvbfhigh1,bins_catvbfhigh1); h0LT_qcd_catvbfhigh_anti->Sumw2();
   TH2F* h0LT_w_catvbfhigh_iso = new TH2F ("h0LT_w_catvbfhigh_iso","h0LT_w_catvbfhigh_iso",binnum_catvbfhigh0,bins_catvbfhigh0,binnum_catvbfhigh1,bins_catvbfhigh1); h0LT_w_catvbfhigh_iso->Sumw2();
   TH2F* h0LT_w_catvbfhigh_anti = new TH2F ("h0LT_w_catvbfhigh_anti","h0LT_w_catvbfhigh_anti",binnum_catvbfhigh0,bins_catvbfhigh0,binnum_catvbfhigh1,bins_catvbfhigh1); h0LT_w_catvbfhigh_anti->Sumw2();
   TH2F* h0J_qcd_catvbfhigh_iso = new TH2F ("h0J_qcd_catvbfhigh_iso","h0J_qcd_catvbfhigh_iso",binnum_catvbfhigh0,bins_catvbfhigh0,binnum_catvbfhigh1,bins_catvbfhigh1); h0J_qcd_catvbfhigh_iso->Sumw2();
   TH2F* h0J_qcd_catvbfhigh_anti = new TH2F ("h0J_qcd_catvbfhigh_anti","h0J_qcd_catvbfhigh_anti",binnum_catvbfhigh0,bins_catvbfhigh0,binnum_catvbfhigh1,bins_catvbfhigh1); h0J_qcd_catvbfhigh_anti->Sumw2();
   TH2F* h0J_w_catvbfhigh_iso = new TH2F ("h0J_w_catvbfhigh_iso","h0J_w_catvbfhigh_iso",binnum_catvbfhigh0,bins_catvbfhigh0,binnum_catvbfhigh1,bins_catvbfhigh1); h0J_w_catvbfhigh_iso->Sumw2();
   TH2F* h0J_w_catvbfhigh_anti = new TH2F ("h0J_w_catvbfhigh_anti","h0J_w_catvbfhigh_anti",binnum_catvbfhigh0,bins_catvbfhigh0,binnum_catvbfhigh1,bins_catvbfhigh1); h0J_w_catvbfhigh_anti->Sumw2();


   string datapath = string(std::getenv("CMSSW_BASE"))+"/src/ComputeFF2018/FFcode/data/";

   reweight::LumiReWeighting* LumiWeights_12;
   LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2018.root").c_str(), (datapath+"pu_distributions_data_2018.root").c_str(), "pileup", "pileup");
   if (year=="2016"){     
     LumiWeights_12 = new reweight::LumiReWeighting((datapath+"MC_Moriond17_PU25ns_V1.root").c_str(), (datapath+"Data_Pileup_2016_271036-284044_80bins.root").c_str(), "pileup", "pileup");
   }
   else if (year=="2017"){
     LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#VBFHToTauTau_M125_13TeV_powheg_pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1#MINIAODSIM", "pileup");
     if (sample=="W1") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3#MINIAODSIM", "pileup");
     else if (sample=="W2") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v5#MINIAODSIM", "pileup");
     else if (sample=="ST_t_antitop") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1#MINIAODSIM", "pileup");
     else if (sample=="DY4") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_v2_94X_mc2017_realistic_v14-v2#MINIAODSIM", "pileup");
     else if (sample=="W") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3#MINIAODSIM", "pileup");
     else if (sample=="DY") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8#RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1#MINIAODSIM", "pileup");
     else if (sample=="WW") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#WW_TuneCP5_13TeV-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1#MINIAODSIM", "pileup");
     else if (sample=="WZ") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#WZ_TuneCP5_13TeV-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1#MINIAODSIM", "pileup");
     else if (sample=="DY4") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1#MINIAODSIM", "pileup");
     else if (sample=="ST_tW_top") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2#MINIAODSIM", "pileup");
     else if (sample=="WminusH125") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#WminusHToTauTau_M125_13TeV_powheg_pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1#MINIAODSIM", "pileup");
     else if (sample=="ST_tW_antitop") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2#MINIAODSIM", "pileup");
     else if (sample=="W3") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1#MINIAODSIM", "pileup");
     else if (sample=="ZH125") LumiWeights_12 = new reweight::LumiReWeighting((datapath+"pu_distributions_mc_2017.root").c_str(), (datapath+"pu_distributions_data_2017.root").c_str(), "pua/#ZHToTauTau_M125_13TeV_powheg_pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1#MINIAODSIM", "pileup");
   }

   TFile fwmc((datapath+"htt_scalefactors_legacy_2018.root").c_str());
   RooWorkspace *wmc = (RooWorkspace*)fwmc.Get("w");
   fwmc.Close();
   if (year=="2017"){
      TFile fwmc((datapath+"htt_scalefactors_legacy_2017.root").c_str());
      wmc = (RooWorkspace*)fwmc.Get("w");
      fwmc.Close();
   }
   if (year=="2016"){
      TFile fwmc((datapath+"htt_scalefactors_legacy_2016.root").c_str());
      wmc = (RooWorkspace*)fwmc.Get("w");
      fwmc.Close();
   }

   TFile frawff("uncorrected_fakefactors_et.root");
   TF1* ff_qcd_0jet=(TF1*) frawff.Get("rawFF_et_qcd_0jet");
   TF1* ff_qcd_1jet=(TF1*) frawff.Get("rawFF_et_qcd_1jet");
   TF1* ff_qcd_2jet=(TF1*) frawff.Get("rawFF_et_qcd_2jet");
   TF1* ff_looseSSqcd_0jet=(TF1*) frawff.Get("rawFF_et_qcd_0jetSSloose");
   TF1* ff_looseSSqcd_1jet=(TF1*) frawff.Get("rawFF_et_qcd_1jetSSloose");
   TF1* ff_looseSSqcd_2jet=(TF1*) frawff.Get("rawFF_et_qcd_2jetSSloose");
   TF1* ff_w_0jet=(TF1*) frawff.Get("rawFF_et_w_0jet");
   TF1* ff_w_1jet=(TF1*) frawff.Get("rawFF_et_w_1jet");
   TF1* ff_w_2jet=(TF1*) frawff.Get("rawFF_et_w_2jet");
   TF1* ff_wmc_0jet=(TF1*) frawff.Get("mc_rawFF_et_w_0jet");
   TF1* ff_wmc_1jet=(TF1*) frawff.Get("mc_rawFF_et_w_1jet");
   TF1* ff_wmc_2jet=(TF1*) frawff.Get("mc_rawFF_et_w_2jet");
   TF1* ff_tt_0jet=(TF1*) frawff.Get("rawFF_et_tt");
   TF1* ff_ttmc_0jet=(TF1*) frawff.Get("mc_rawFF_et_tt");

   TFile fmvisclosure ("FF_corrections_1.root");
   TF1* mvisclosure_wmc=(TF1*) fmvisclosure.Get("closure_mvis_et_wmc");
   TF1* mvisclosure_qcd=(TF1*) fmvisclosure.Get("closure_mvis_et_qcd");
   TF1* mvisclosure_w=(TF1*) fmvisclosure.Get("closure_mvis_et_w");
   TF1* mvisclosure_tt=(TF1*) fmvisclosure.Get("closure_mvis_et_ttmc");
   TF1* mvisclosure_qcd0=(TF1*) fmvisclosure.Get("closure_mvis_et_0jet_qcd");
   TF1* mvisclosure_w0=(TF1*) fmvisclosure.Get("closure_mvis_et_0jet_w");
   TF1* mvisclosure_qcd1=(TF1*) fmvisclosure.Get("closure_mvis_et_1jet_qcd");
   TF1* mvisclosure_w1=(TF1*) fmvisclosure.Get("closure_mvis_et_1jet_w");
   TF1* mvisclosure_qcd2=(TF1*) fmvisclosure.Get("closure_mvis_et_2jet_qcd");
   TF1* mvisclosure_w2=(TF1*) fmvisclosure.Get("closure_mvis_et_2jet_w");

   ScaleFactor * myScaleFactor_trgEle3235 = new ScaleFactor();
   myScaleFactor_trgEle3235->init_ScaleFactor(string(std::getenv("CMSSW_BASE"))+"/src/LeptonEfficiencies/Electron/Run2018/Electron_Run2018_Ele32orEle35.root");
   ScaleFactor * myScaleFactor_trgEle24 = new ScaleFactor();
   myScaleFactor_trgEle24->init_ScaleFactor(string(std::getenv("CMSSW_BASE"))+"/src/LeptonEfficiencies/Electron/Run2018/Electron_Run2018_Ele24.root");
   ScaleFactor * myScaleFactor_IdIso = new ScaleFactor();
   myScaleFactor_IdIso->init_ScaleFactor(string(std::getenv("CMSSW_BASE"))+"/src/LeptonEfficiencies/Electron/Run2018/Electron_Run2018_IdIso.root");

   if (year=="2017"){
     myScaleFactor_trgEle3235->init_ScaleFactor(string(std::getenv("CMSSW_BASE"))+"/src/LeptonEfficiencies/Electron/Run2017/Electron_Ele32orEle35.root");
     myScaleFactor_trgEle24->init_ScaleFactor(string(std::getenv("CMSSW_BASE"))+"/src/LeptonEfficiencies/Electron/Run2017/Electron_EleTau_Ele24.root");
     myScaleFactor_IdIso->init_ScaleFactor(string(std::getenv("CMSSW_BASE"))+"/src/LeptonEfficiencies/Electron/Run2017/Electron_IdIso_IsoLt0.15_IsoID_eff.root");
   }

   TauTriggerSFs2017* etsf=new TauTriggerSFs2017(string(std::getenv("CMSSW_BASE"))+"/src/TauAnalysisTools/TauTriggerSFs/data/tauTriggerEfficiencies2018.root","etau", "2018", "tight", "MVAv2");
   if (year=="2017") etsf=new TauTriggerSFs2017(string(std::getenv("CMSSW_BASE"))+"/src/TauAnalysisTools/TauTriggerSFs/data/tauTriggerEfficiencies2017.root","etau", "2017", "tight", "MVAv2");

   TFile *fZ=new TFile((datapath+"zpt_weights_2016_BtoH.root").c_str());
   TH2F *histZ=(TH2F*) fZ->Get("zptmass_histo");
   ScaleFactor * myScaleFactor_trgEle25 = new ScaleFactor();
   myScaleFactor_trgEle25->init_ScaleFactor(string(std::getenv("CMSSW_BASE"))+"/src/LeptonEfficiencies/Electron/Run2016_legacy/Electron_Run2016_legacy_Ele25.root");
   if (year=="2016") myScaleFactor_IdIso->init_ScaleFactor(string(std::getenv("CMSSW_BASE"))+"/src/LeptonEfficiencies/Electron/Run2016_legacy/Electron_Run2016_legacy_IdIso.root");

   TauIDSFTool * theSFTool;
   if (year == "2016") theSFTool = new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSjet","Medium");
   else if (year == "2017") theSFTool = new TauIDSFTool("2017ReReco","DeepTau2017v2p1VSjet","Medium");
   else theSFTool = new TauIDSFTool("2018ReReco","DeepTau2017v2p1VSjet", "Medium");

   Int_t nentries_wtn = (Int_t) arbre->GetEntries();
   for (Int_t i = 0; i < nentries_wtn; i++) {
        arbre->GetEntry(i);	
        if (i % 10000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
        fflush(stdout);

	if (fabs(eta_1)>2.1) continue;
        if (fabs(eta_2)>2.3) continue;

        if (Flag_goodVertices) continue;
        if (Flag_globalSuperTightHalo2016Filter) continue;
        if (Flag_HBHENoiseFilter) continue;
        if (Flag_HBHENoiseIsoFilter) continue;
        if (Flag_EcalDeadCellTriggerPrimitiveFilter) continue;
        if (Flag_BadPFMuonFilter) continue;
        if ((sample=="data_obs" or sample=="embedded") && Flag_eeBadScFilter) continue;
        if (Flag_ecalBadCalibReducedMINIAODFilter) continue;

        bool trigger35=(passEle35 && pt_1>33 && matchEle35_1 && filterEle35_1);
        bool trigger32=(passEle32 && pt_1>33 && matchEle32_1 && filterEle32_1);
        bool trigger2430=(passEle24Tau30 && matchEle24Tau30_1 && filterEle24Tau30_1 && matchEle24Tau30_2 && filterEle24Tau30_2 && pt_1>25 && pt_2>35 && fabs(eta_2)<2.1 && pt_1<=33);
	bool trigger27=false;
	bool trigger25=false;
        bool trigger2430HPS=(passEle24HPSTau30 && matchEle24HPSTau30_1 && filterEle24HPSTau30_1 && matchEle24HPSTau30_2 && filterEle24HPSTau30_2 && pt_1>25 && pt_2>35 && fabs(eta_2)<2.1 && pt_1<=33);
	if (sample=="embedded"){
	   if (fabs(eta_1)>=1.479){
		trigger35=(pt_1>33); trigger32=(pt_1>33); trigger2430HPS=(fabs(eta_2)<2.1 && pt_1>25 && pt_2>35 && pt_1<=33);
	   }
	   if (fabs(eta_1)<1.479){   
                trigger35=(pt_1>33 && matchEmbFilter_Ele35_1); trigger32=(pt_1>33 && (matchEmbFilter_Ele32DoubleL1v1_1 or matchEmbFilter_Ele32DoubleL1v2_1 or matchEmbFilter_Ele32_1)); trigger2430HPS=(matchEmbFilter_Ele24Tau30_1 && matchEmbFilter_Ele24Tau30_2 && fabs(eta_2)<2.1 && pt_1>25 && pt_2>35 && pt_1<=33);
           }
	}
        if (year=="2018" && sample=="data_obs" && run<317509 && !trigger2430 && !trigger32 && !trigger35) continue;
        if (year=="2018" && sample=="data_obs" && run>=317509 && !trigger2430HPS && !trigger32 && !trigger35) continue;
        if (year=="2018" && sample!="data_obs" && !trigger32 && !trigger35 && !trigger2430HPS) continue;

        if (year=="2017"){
           trigger35=(passEle35 && pt_1>28 && matchEle35_1 && filterEle35_1);
           trigger32=(passEle32 && pt_1>28 && matchEle32_1 && filterEle32_1);
           trigger27=(passEle27 && pt_1>28 && matchEle27_1 && filterEle27_1);
           trigger2430=(passEle24Tau30 && matchEle24Tau30_1 && filterEle24Tau30_1 && matchEle24Tau30_2 && filterEle24Tau30_2 && pt_1>25 && pt_2>35 && pt_1<28 && fabs(eta_2)<2.1);
           if (sample=="embedded" && fabs(eta_1)<1.479) trigger2430=(passEle24Tau30 && pt_1>25 && pt_2>35 && pt_1<28 && fabs(eta_2)<2.1);
           if (sample=="embedded" && fabs(eta_1)>1.479 && (pt_1>28 or fabs(eta_2)<2.1)){ trigger35=true; trigger32=true; trigger2430=true;}
           if (sample=="embedded" && fabs(eta_1)>1.479 && pt_1<28 && (pt_2<35 or fabs(eta_2)>2.1)){ trigger35=false; trigger32=false; trigger2430=false;}
           if (!trigger27 && !trigger32 && !trigger35 && !trigger2430) continue;
        }
        else if (year=="2016"){
           trigger25=(passEle25 && matchEle25_1 && filterEle25_1 && pt_1>26);
           if (!trigger25) continue;
        }
	


        // Deep Medium
        if (!byTightDeepVSe_2 or !byVLooseDeepVSmu_2) continue;
        float signalRegion=(byMediumDeepVSjet_2);
        float antiisoRegion=(byVVVLooseDeepVSjet_2 && !byMediumDeepVSjet_2);

	TLorentzVector mytau; 
	mytau.SetPtEtaPhiM(pt_2,eta_2,phi_2,m_2);
        TLorentzVector myele;
        myele.SetPtEtaPhiM(pt_1,eta_1,phi_1,m_1);

	if (myele.DeltaR(mytau)<0.5) continue;

	if (year=="2018"){
           if (sample=="W"){
               weight=51.749;
               if (numGenJets==1) weight=10.170;
               else if (numGenJets==2) weight=4.51855;
               else if (numGenJets==3) weight=3.07747;
               else if (numGenJets==4) weight=3.2113;
           }

           if (sample=="DY"){
               weight=3.6234;
               if (numGenJets==1) weight=0.6298;
               else if (numGenJets==2) weight=0.5521;
               else if (numGenJets==3) weight=0.5995;
               else if (numGenJets==4) weight=0.8211;
           }
	}
        else if (year=="2017"){
           if (sample=="W"){
               weight=23.8336;
               if (numGenJets==1) weight=3.1468;
               else if (numGenJets==2) weight=4.087;
               else if (numGenJets==3) weight=2.253;
               else if (numGenJets==4) weight=2.1954;
           }
           if (sample=="DY"){
               weight=2.61364;
               if (numGenJets==1) weight=0.71866;
               else if (numGenJets==2) weight=0.9380;
               else if (numGenJets==3) weight=1.6643;
               else if (numGenJets==4) weight=0.2359;
           }
        }
        else if (year=="2016"){
           if (sample=="W"){
               weight=25.3918;
               if (numGenJets==1) weight=5.76634;
               else if (numGenJets==2) weight=1.7906;
               else if (numGenJets==3) weight=0.67907;
               else if (numGenJets==4) weight=1.84847;
           }
           if (sample=="DY"){
               weight=1.49237;
               if (numGenJets==1) weight=0.47595;
               else if (numGenJets==2) weight=0.4931;
               else if (numGenJets==3) weight=0.5113;
               else if (numGenJets==4) weight=0.69041;
           }
        }

        bool is_includedInEmbedded=false;
        //if ((name.find("125")>100 && sample!="data_obs" && sample!="embedded") && gen_match_1>2 && gen_match_1<6 && gen_match_2>2 && gen_match_2<6) is_includedInEmbedded=true; // remove overlap with embedded samples
        bool isT=(!is_includedInEmbedded && gen_match_2==5);
        bool isL=(!is_includedInEmbedded && gen_match_2<5);

        float aweight=genweight*weight*LumiWeights_12->weight(npu);
        if (sample=="embedded") aweight=genweight;
        if (sample!="embedded" && byMediumDeepVSjet_2 && sample!="data_obs" && gen_match_2==5) aweight=aweight*theSFTool->getSFvsPT(mytau.Pt()); // deep M

        //if (byLooseDeepVSjet_2 && sample!="embedded" && sample!="data_obs" && gen_match_2==5) aweight=aweight*0.86;
        if (sample=="embedded") aweight=aweight*0.97;
	/*if (gen_match_2==2 or gen_match_2==4){
	   if (fabs(mytau.Eta())<0.4) aweight=aweight*1.06;
           else if (fabs(mytau.Eta())<0.8) aweight=aweight*1.02;
           else if (fabs(mytau.Eta())<1.2) aweight=aweight*1.10;
           else if (fabs(mytau.Eta())<1.7) aweight=aweight*1.03;
           else if (fabs(mytau.Eta())<2.3) aweight=aweight*1.94;
	}
        if (gen_match_2==1 or gen_match_2==3){
           if (fabs(mytau.Eta())<1.460) aweight=aweight*1.80;
           else if (fabs(mytau.Eta())>1.558) aweight=aweight*1.53;
	}*/

	TLorentzVector mymet;
	mymet.SetPtEtaPhiM(met,0,metphi,0);

	if (sample=="data_obs") aweight=1.0;


	// Top pT reweighting
        float topfactor=1.0;
        if (name=="TT" or name=="TTMC"){
           float pttop1=pt_top1;
           if (pttop1>400) pttop1=400;
           float pttop2=pt_top2;
           if (pttop2>400) pttop2=400;
           topfactor=sqrt(exp(0.0615-0.0005*pttop1)*exp(0.0615-0.0005*pttop2));
           aweight*=topfactor;
        }

        float zptweight=1.0;
	if (year=="2018" && sample!="embedded" && sample!="data_obs"){
          wmc->var("z_gen_mass")->setVal(genM);
          wmc->var("z_gen_pt")->setVal(genpT);
	  zptweight=wmc->function("zptmass_weight_nom")->getVal();
	  if (sample=="DY") aweight=aweight*zptweight;
	  aweight=aweight*myScaleFactor_IdIso->get_ScaleFactor(pt_1,eta_1);
             int mydm=l2_decayMode;
             if (mydm>10) mydm=10;
	  if (trigger32 or trigger35) aweight=aweight*myScaleFactor_trgEle3235->get_ScaleFactor(pt_1,eta_1);
	  else aweight=aweight*myScaleFactor_trgEle24->get_ScaleFactor(pt_1,eta_1)*etsf->getTriggerScaleFactor(mytau.Pt(), mytau.Eta(), mytau.Phi(), mydm);
          aweight=aweight*bweight;
	}	

        if (year=="2017" && sample!="embedded" && sample!="data_obs"){
          wmc->var("e_pt")->setVal(myele.Pt());
          wmc->var("e_eta")->setVal(myele.Eta());
          wmc->var("e_iso")->setVal(iso_1);
          wmc->var("z_gen_mass")->setVal(genM);
          wmc->var("z_gen_pt")->setVal(genpT);
          aweight=aweight*wmc->function("e_trk_ratio")->getVal();
          aweight=aweight*myScaleFactor_IdIso->get_ScaleFactor(pt_1,eta_1);
          aweight=aweight*0.991;//z vtx electron HLT
          zptweight=wmc->function("zptmass_weight_nom")->getVal();
          if (sample=="DY") aweight=aweight*zptweight;
          int dm=l2_decayMode;
          if (dm>10) dm=10;
          if (trigger27 or trigger32 or trigger35) aweight=aweight*wmc->function("e_trg27_trg32_trg35_kit_data")->getVal()/wmc->function("e_trg27_trg32_trg35_kit_mc")->getVal();
          else aweight=aweight*(wmc->function("e_trg_EleTau_Ele24Leg_desy_data")->getVal()/wmc->function("e_trg_EleTau_Ele24Leg_desy_mc")->getVal())*etsf->getTriggerScaleFactor(mytau.Pt(), mytau.Eta(), mytau.Phi(), dm);
          aweight=aweight*bweight*prefiring_weight;
        }

        if (year=="2016" && sample!="embedded" && sample!="data_obs"){
          aweight=aweight*myScaleFactor_IdIso->get_ScaleFactor(pt_1,eta_1);
          aweight=aweight*myScaleFactor_trgEle25->get_ScaleFactor(pt_1,eta_1);
          wmc->var("z_gen_mass")->setVal(genM);
          wmc->var("z_gen_pt")->setVal(genpT);
          zptweight=wmc->function("zptmass_weight_nom")->getVal();
          if (sample=="DY") aweight=aweight*zptweight;
          wmc->var("e_pt")->setVal(myele.Pt());
          wmc->var("e_eta")->setVal(myele.Eta());
          wmc->var("e_iso")->setVal(iso_1);
          aweight=aweight*wmc->function("e_trk_ratio")->getVal();
          aweight=aweight*bweight*prefiring_weight;
        }

	float mt=TMass_F(myele.Pt(),mymet.Pt(),myele.Px(),mymet.Px(),myele.Py(),mymet.Py());

	//************************* Fill histograms **********************
           if (mytau.Pt()<30) continue;
	   float weight2=1.0;
	   float myvar=(myele+mytau).M();
	   //if (myvar>300) myvar=299;

	  float ff_qcd=get_raw_FF(mytau.Pt(),ff_qcd_0jet);
	  if (njets==1) ff_qcd=get_raw_FF(mytau.Pt(),ff_qcd_1jet);
          else if (njets>1) ff_qcd=get_raw_FF(mytau.Pt(),ff_qcd_2jet);
          float ff_looseSSqcd=get_raw_FF(mytau.Pt(),ff_looseSSqcd_0jet);
          if (njets==1) ff_looseSSqcd=get_raw_FF(mytau.Pt(),ff_looseSSqcd_1jet);
          else if (njets>1) ff_looseSSqcd=get_raw_FF(mytau.Pt(),ff_looseSSqcd_2jet);
          float ff_w=get_raw_FF(mytau.Pt(),ff_w_0jet);
          if (njets==1) ff_w=get_raw_FF(mytau.Pt(),ff_w_1jet);
          else if (njets>1) ff_w=get_raw_FF(mytau.Pt(),ff_w_2jet);
          float ff_tt=get_raw_FF(mytau.Pt(),ff_tt_0jet);

          if (apply_corrections=="yes"){
             ff_qcd=get_raw_FF(mytau.Pt(),ff_qcd_0jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_qcd0);
             if (njets==1) ff_qcd=get_raw_FF(mytau.Pt(),ff_qcd_1jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_qcd1);
             else if (njets>1) ff_qcd=get_raw_FF(mytau.Pt(),ff_qcd_2jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_qcd2);
             ff_w=get_raw_FF(mytau.Pt(),ff_w_0jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_w0);
             if (njets==1) ff_w=get_raw_FF(mytau.Pt(),ff_w_1jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_w1);
             else if (njets>1) ff_w=get_raw_FF(mytau.Pt(),ff_w_2jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_w2);
             ff_tt=get_raw_FF(mytau.Pt(),ff_tt_0jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_tt);
	  }

	  if (name=="WMC"){
	     ff_w=get_raw_FF(mytau.Pt(),ff_wmc_0jet);
             if (njets==1) ff_w=get_raw_FF(mytau.Pt(),ff_wmc_1jet);
             else if (njets>1) ff_w=get_raw_FF(mytau.Pt(),ff_wmc_2jet);
	     mt=100;
	  }

          if (name=="TTMC"){
             ff_tt=get_raw_FF(mytau.Pt(),ff_ttmc_0jet);
          }	  

          if (name=="WMC2"){	    	    
             ff_w=get_raw_FF(mytau.Pt(),ff_wmc_0jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_wmc);
             if (njets==1) ff_w=get_raw_FF(mytau.Pt(),ff_wmc_1jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_wmc);
             else if (njets>1) ff_w=get_raw_FF(mytau.Pt(),ff_wmc_2jet)*get_mvis_closure((myele+mytau).M(),mvisclosure_wmc);
          }	  

           if (!is_includedInEmbedded){

	     // W MC
	     if (signalRegion && iso_1<0.15 && nbtag==0 && q_1*q_2<0)
                  hmt_w_iso->Fill(mt,aweight*weight2);
             if (antiisoRegion && iso_1<0.15 && nbtag==0 && q_1*q_2<0)
                  hmt_w_anti->Fill(mt,aweight*weight2*ff_w);

	     if (isL or isT){

	       // QCD
               if (signalRegion && iso_1>0.02 && iso_1<0.15 && nbtag==0 && mt<50 && q_1*q_2>0){
                  if (njets==0) h0LT_qcd_mvis_iso->Fill(myvar,aweight*weight2);
                  else if (njets==1) h1LT_qcd_mvis_iso->Fill(myvar,aweight*weight2);
                  else if (njets>1) h2LT_qcd_mvis_iso->Fill(myvar,aweight*weight2);
                  h0LT_qcd_iso->Fill(myvar,aweight*weight2);
                  h0LT_qcd_taupt_iso->Fill(mytau.Pt(),aweight*weight2);
                  if (njets>0) h0LT_qcd_pth_iso->Fill((myele+mymet+mytau).Pt(),aweight*weight2);
                  if (njets>1) h0LT_qcd_mjj_iso->Fill(mjj,aweight*weight2);
                  if (njets==0 && (myele+mytau+mymet).Pt()<10) h0LT_qcd_cat0jetlow_iso->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2);
                  if (njets==0 && (myele+mytau+mymet).Pt()>10) h0LT_qcd_cat0jethigh_iso->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2);
                  if (njets==1) h0LT_qcd_catboosted1_iso->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2);
                  if (njets>1 && mjj<350) h0LT_qcd_catboosted2_iso->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()<200) h0LT_qcd_catvbflow_iso->Fill(mjj,m_sv,aweight*weight2);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()>200) h0LT_qcd_catvbfhigh_iso->Fill(mjj,m_sv,aweight*weight2);
               }
               if (antiisoRegion && iso_1>0.02 && iso_1<0.15 && nbtag==0 && mt<50 && q_1*q_2>0){
                  if (njets==0) h0LT_qcd_mvis_anti->Fill(myvar,aweight*weight2*ff_qcd);
                  else if (njets==1) h1LT_qcd_mvis_anti->Fill(myvar,aweight*weight2*ff_qcd);
                  else if (njets>1) h2LT_qcd_mvis_anti->Fill(myvar,aweight*weight2*ff_qcd);
                  h0LT_qcd_anti->Fill(myvar,aweight*weight2*ff_qcd);
                  h0LT_qcd_taupt_anti->Fill(mytau.Pt(),aweight*weight2*ff_qcd);
                  if (njets>0) h0LT_qcd_pth_anti->Fill((myele+mymet+mytau).Pt(),aweight*weight2*ff_qcd);
                  if (njets>1) h0LT_qcd_mjj_anti->Fill(mjj,aweight*weight2*ff_qcd);
                  if (njets==0 && (myele+mytau+mymet).Pt()<10) h0LT_qcd_cat0jetlow_anti->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2*ff_qcd);
                  if (njets==0 && (myele+mytau+mymet).Pt()>10) h0LT_qcd_cat0jethigh_anti->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2*ff_qcd);
                  if (njets==1) h0LT_qcd_catboosted1_anti->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2*ff_qcd);
                  if (njets>1 && mjj<350) h0LT_qcd_catboosted2_anti->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2*ff_qcd);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()<200) h0LT_qcd_catvbflow_anti->Fill(mjj,m_sv,aweight*weight2*ff_qcd);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()>200) h0LT_qcd_catvbfhigh_anti->Fill(mjj,m_sv,aweight*weight2*ff_qcd);
               }

	       // QCD loose
               if (signalRegion && iso_1>0.15 && iso_1<0.25 && nbtag==0 && mt<50 && q_1*q_2>0)
                  h0SSlooseLT_qcd_iso->Fill(myvar,aweight*weight2);
               if (antiisoRegion && iso_1>0.15 && iso_1<0.25 && nbtag==0 && mt<50 && q_1*q_2>0)
                  h0SSlooseLT_qcd_anti->Fill(myvar,aweight*weight2*ff_looseSSqcd);


	       // W
               if (signalRegion && iso_1<0.15 && nbtag==0 && mt>70 && q_1*q_2<0){
		  if (njets==0) h0LT_w_mvis_iso->Fill(myvar,aweight*weight2);
                  else if (njets==1) h1LT_w_mvis_iso->Fill(myvar,aweight*weight2);
                  else if (njets>1) h2LT_w_mvis_iso->Fill(myvar,aweight*weight2);
                  h0LT_w_iso->Fill(myvar,aweight*weight2);
                  h0LT_w_taupt_iso->Fill(mytau.Pt(),aweight*weight2);
                  if (njets>0) h0LT_w_pth_iso->Fill((myele+mymet+mytau).Pt(),aweight*weight2);
                  if (njets>1) h0LT_w_mjj_iso->Fill(mjj,aweight*weight2);
		  if (njets==0 && (myele+mytau+mymet).Pt()<10) h0LT_w_cat0jetlow_iso->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2);
		  if (njets==0 && (myele+mytau+mymet).Pt()>10) h0LT_w_cat0jethigh_iso->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2);
		  if (njets==1) h0LT_w_catboosted1_iso->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2);
		  if (njets>1 && mjj<350) h0LT_w_catboosted2_iso->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2);
		  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()<200) h0LT_w_catvbflow_iso->Fill(mjj,m_sv,aweight*weight2);
		  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()>200) h0LT_w_catvbfhigh_iso->Fill(mjj,m_sv,aweight*weight2);
	       }
               if (antiisoRegion && iso_1<0.15 && nbtag==0 && mt>70 && q_1*q_2<0){
		  if (njets==0) h0LT_w_mvis_anti->Fill(myvar,aweight*weight2*ff_w);
                  else if (njets==1) h1LT_w_mvis_anti->Fill(myvar,aweight*weight2*ff_w);
                  else if (njets>1) h2LT_w_mvis_anti->Fill(myvar,aweight*weight2*ff_w);
                  h0LT_w_anti->Fill(myvar,aweight*weight2*ff_w);
                  h0LT_w_taupt_anti->Fill(mytau.Pt(),aweight*weight2*ff_w);
                  if (njets>0) h0LT_w_pth_anti->Fill((myele+mymet+mytau).Pt(),aweight*weight2*ff_w);
                  if (njets>1) h0LT_w_mjj_anti->Fill(mjj,aweight*weight2*ff_w);
                  if (njets==0 && (myele+mytau+mymet).Pt()<10) h0LT_w_cat0jetlow_anti->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2*ff_w);
                  if (njets==0 && (myele+mytau+mymet).Pt()>10) h0LT_w_cat0jethigh_anti->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2*ff_w);
                  if (njets==1) h0LT_w_catboosted1_anti->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2*ff_w);
                  if (njets>1 && mjj<350) h0LT_w_catboosted2_anti->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2*ff_w);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()<200) h0LT_w_catvbflow_anti->Fill(mjj,m_sv,aweight*weight2*ff_w);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()>200) h0LT_w_catvbfhigh_anti->Fill(mjj,m_sv,aweight*weight2*ff_w);
	       }

	       // TT
               if (signalRegion && iso_1<0.15 && nbtag>0 && mt<50 && q_1*q_2<0)
                  h0LT_tt_iso->Fill(myvar,aweight*weight2);
               if (antiisoRegion && iso_1<0.15 && nbtag>0 && mt<50 && q_1*q_2<0)
                  h0LT_tt_anti->Fill(myvar,aweight*weight2*ff_tt);

	    }
	   else{
               if (signalRegion && iso_1>0.02 && iso_1<0.15 && nbtag==0 && mt<50 && q_1*q_2>0){
                  if (njets==0) h0J_qcd_mvis_iso->Fill(myvar,aweight*weight2);
                  else if (njets==1) h1J_qcd_mvis_iso->Fill(myvar,aweight*weight2);
                  else if (njets>1) h2J_qcd_mvis_iso->Fill(myvar,aweight*weight2);
                  h0J_qcd_iso->Fill(myvar,aweight*weight2);
                  h0J_qcd_taupt_iso->Fill(mytau.Pt(),aweight*weight2);
                  if (njets>0) h0J_qcd_pth_iso->Fill((myele+mymet+mytau).Pt(),aweight*weight2);
                  if (njets>1) h0J_qcd_mjj_iso->Fill(mjj,aweight*weight2);
                  if (njets==0 && (myele+mytau+mymet).Pt()<10) h0J_qcd_cat0jetlow_iso->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2);
                  if (njets==0 && (myele+mytau+mymet).Pt()>200) h0J_qcd_cat0jethigh_iso->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2);
                  if (njets==1) h0J_qcd_catboosted1_iso->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2);
                  if (njets>1 && mjj<350) h0J_qcd_catboosted2_iso->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()<200) h0J_qcd_catvbflow_iso->Fill(mjj,m_sv,aweight*weight2);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()>200) h0J_qcd_catvbfhigh_iso->Fill(mjj,m_sv,aweight*weight2);

               }
               if (antiisoRegion && iso_1>0.02 && iso_1<0.15 && nbtag==0 && mt<50 && q_1*q_2>0){
                  if (njets==0) h0J_qcd_mvis_anti->Fill(myvar,aweight*weight2*ff_qcd);
                  else if (njets==1) h1J_qcd_mvis_anti->Fill(myvar,aweight*weight2*ff_qcd);
                  else if (njets>1) h2J_qcd_mvis_anti->Fill(myvar,aweight*weight2*ff_qcd);
                  h0J_qcd_anti->Fill(myvar,aweight*weight2*ff_qcd);
                  h0J_qcd_taupt_anti->Fill(mytau.Pt(),aweight*weight2*ff_qcd);
                  if (njets>0) h0J_qcd_pth_anti->Fill((myele+mymet+mytau).Pt(),aweight*weight2*ff_qcd);
                  if (njets>1) h0J_qcd_mjj_anti->Fill(mjj,aweight*weight2*ff_qcd);
                  if (njets==0 && (myele+mytau+mymet).Pt()<10) h0J_qcd_cat0jetlow_anti->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2*ff_qcd);
                  if (njets==0 && (myele+mytau+mymet).Pt()>200) h0J_qcd_cat0jethigh_anti->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2*ff_qcd);
                  if (njets==1) h0J_qcd_catboosted1_anti->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2*ff_qcd);
                  if (njets>1 && mjj<350) h0J_qcd_catboosted2_anti->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2*ff_qcd);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()<200) h0J_qcd_catvbflow_anti->Fill(mjj,m_sv,aweight*weight2*ff_qcd);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()>200) h0J_qcd_catvbfhigh_anti->Fill(mjj,m_sv,aweight*weight2*ff_qcd);
	       }

	       // QCD loose
               if (signalRegion && iso_1>0.15 && iso_1<0.25 && nbtag==0 && mt<50 && q_1*q_2>0)
                  h0SSlooseJ_qcd_iso->Fill(myvar,aweight*weight2);
               if (antiisoRegion && iso_1>0.15 && iso_1<0.25 && nbtag==0 && mt<50 && q_1*q_2>0)
                  h0SSlooseJ_qcd_anti->Fill(myvar,aweight*weight2*ff_looseSSqcd);

	       // W 
               if (signalRegion && iso_1<0.15 && nbtag==0 && mt>70 && q_1*q_2<0){
		  if (njets==0) h0J_w_mvis_iso->Fill(myvar,aweight*weight2);
                  else if (njets==1) h1J_w_mvis_iso->Fill(myvar,aweight*weight2);
                  else if (njets>1) h2J_w_mvis_iso->Fill(myvar,aweight*weight2);
                  h0J_w_iso->Fill(myvar,aweight*weight2);
                  h0J_w_taupt_iso->Fill(mytau.Pt(),aweight*weight2);
                  if (njets>0) h0J_w_pth_iso->Fill((myele+mymet+mytau).Pt(),aweight*weight2);
                  if (njets>1) h0J_w_mjj_iso->Fill(mjj,aweight*weight2);
		  if (njets==0 && (myele+mytau+mymet).Pt()<10) h0J_w_cat0jetlow_iso->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2);
		  if (njets==0 && (myele+mytau+mymet).Pt()>200) h0J_w_cat0jethigh_iso->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2);
		  if (njets==1) h0J_w_catboosted1_iso->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2);
		  if (njets>1 && mjj<350) h0J_w_catboosted2_iso->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2);
		  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()<200) h0J_w_catvbflow_iso->Fill(mjj,m_sv,aweight*weight2);
		  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()>200) h0J_w_catvbfhigh_iso->Fill(mjj,m_sv,aweight*weight2);

	       }
               if (antiisoRegion && iso_1<0.15 && nbtag==0 && mt>70 && q_1*q_2<0){
		  if (njets==0) h0J_w_mvis_anti->Fill(myvar,aweight*weight2*ff_w);
                  else if (njets==1) h1J_w_mvis_anti->Fill(myvar,aweight*weight2*ff_w);
                  else if (njets>1) h2J_w_mvis_anti->Fill(myvar,aweight*weight2*ff_w);
                  h0J_w_anti->Fill(myvar,aweight*weight2*ff_w);
                  h0J_w_taupt_anti->Fill(mytau.Pt(),aweight*weight2*ff_w);
                  if (njets>0) h0J_w_pth_anti->Fill((myele+mymet+mytau).Pt(),aweight*weight2*ff_w);
                  if (njets>1) h0J_w_mjj_anti->Fill(mjj,aweight*weight2*ff_w);
                  if (njets==0 && (myele+mytau+mymet).Pt()<10) h0J_w_cat0jetlow_anti->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2*ff_w);
                  if (njets==0 && (myele+mytau+mymet).Pt()>200) h0J_w_cat0jethigh_anti->Fill(mytau.Pt(),(myele+mytau).M(),aweight*weight2*ff_w);
                  if (njets==1) h0J_w_catboosted1_anti->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2*ff_w);
                  if (njets>1 && mjj<350) h0J_w_catboosted2_anti->Fill((myele+mytau+mymet).Pt(),m_sv,aweight*weight2*ff_w);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()<200) h0J_w_catvbflow_anti->Fill(mjj,m_sv,aweight*weight2*ff_w);
                  if (njets>1 && mjj>350 && (myele+mytau+mymet).Pt()>200) h0J_w_catvbfhigh_anti->Fill(mjj,m_sv,aweight*weight2*ff_w);

	       }

	       // TT
               if (signalRegion && iso_1<0.15 && nbtag>0 && mt<50 && q_1*q_2<0)
                  h0J_tt_iso->Fill(myvar,aweight*weight2);
               if (antiisoRegion && iso_1<0.15 && nbtag>0 && mt<50 && q_1*q_2<0)
                  h0J_tt_anti->Fill(myvar,aweight*weight2*ff_tt);
            }

           }

    } // end of loop over events
    TFile *fout = TFile::Open(output.c_str(), "RECREATE");
    fout->cd();

    TString postfixJ="J";
    TString postfixLT="LT";

    TDirectory *d0_qcd_mvis_iso =fout->mkdir("et_0jet_qcd_mvis_iso");
    d0_qcd_mvis_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_mvis_iso->SetName(name.c_str());
      h0LT_qcd_mvis_iso->Add(h0J_qcd_mvis_iso);
      h0LT_qcd_mvis_iso->Write();
    }
    else{
      h0LT_qcd_mvis_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_mvis_iso->Write();
      h0J_qcd_mvis_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_mvis_iso->Write();
    }
    TDirectory *d0_qcd_mvis_anti =fout->mkdir("et_0jet_qcd_mvis_anti");
    d0_qcd_mvis_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_mvis_anti->SetName(name.c_str());
      h0LT_qcd_mvis_anti->Add(h0J_qcd_mvis_anti);
      h0LT_qcd_mvis_anti->Write();
    }
    else{
      h0LT_qcd_mvis_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_mvis_anti->Write();
      h0J_qcd_mvis_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_mvis_anti->Write();
    }

    TDirectory *d1_qcd_mvis_iso =fout->mkdir("et_1jet_qcd_mvis_iso");
    d1_qcd_mvis_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h1LT_qcd_mvis_iso->SetName(name.c_str());
      h1LT_qcd_mvis_iso->Add(h1J_qcd_mvis_iso);
      h1LT_qcd_mvis_iso->Write();
    }
    else{
      h1LT_qcd_mvis_iso->SetName(name.c_str()+postfixLT);
      h1LT_qcd_mvis_iso->Write();
      h1J_qcd_mvis_iso->SetName(name.c_str()+postfixJ);
      h1J_qcd_mvis_iso->Write();
    }
    TDirectory *d1_qcd_mvis_anti =fout->mkdir("et_1jet_qcd_mvis_anti");
    d1_qcd_mvis_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h1LT_qcd_mvis_anti->SetName(name.c_str());
      h1LT_qcd_mvis_anti->Add(h1J_qcd_mvis_anti);
      h1LT_qcd_mvis_anti->Write();
    }
    else{
      h1LT_qcd_mvis_anti->SetName(name.c_str()+postfixLT);
      h1LT_qcd_mvis_anti->Write();
      h1J_qcd_mvis_anti->SetName(name.c_str()+postfixJ);
      h1J_qcd_mvis_anti->Write();
    }

    TDirectory *d2_qcd_mvis_iso =fout->mkdir("et_2jet_qcd_mvis_iso");
    d2_qcd_mvis_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h2LT_qcd_mvis_iso->SetName(name.c_str());
      h2LT_qcd_mvis_iso->Add(h2J_qcd_mvis_iso);
      h2LT_qcd_mvis_iso->Write();
    }
    else{
      h2LT_qcd_mvis_iso->SetName(name.c_str()+postfixLT);
      h2LT_qcd_mvis_iso->Write();
      h2J_qcd_mvis_iso->SetName(name.c_str()+postfixJ);
      h2J_qcd_mvis_iso->Write();
    }
    TDirectory *d2_qcd_mvis_anti =fout->mkdir("et_2jet_qcd_mvis_anti");
    d2_qcd_mvis_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h2LT_qcd_mvis_anti->SetName(name.c_str());
      h2LT_qcd_mvis_anti->Add(h2J_qcd_mvis_anti);
      h2LT_qcd_mvis_anti->Write();
    }
    else{
      h2LT_qcd_mvis_anti->SetName(name.c_str()+postfixLT);
      h2LT_qcd_mvis_anti->Write();
      h2J_qcd_mvis_anti->SetName(name.c_str()+postfixJ);
      h2J_qcd_mvis_anti->Write();
    }

    TDirectory *d0_qcd_iso =fout->mkdir("et_0jet_qcd_iso");
    d0_qcd_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_iso->SetName(name.c_str());
      h0LT_qcd_iso->Add(h0J_qcd_iso);
      h0LT_qcd_iso->Write();
    }
    else{
      h0LT_qcd_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_iso->Write();
      h0J_qcd_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_iso->Write();
    }
    TDirectory *d0_qcd_anti =fout->mkdir("et_0jet_qcd_anti");
    d0_qcd_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_anti->SetName(name.c_str());
      h0LT_qcd_anti->Add(h0J_qcd_anti);
      h0LT_qcd_anti->Write();
    }
    else{
      h0LT_qcd_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_anti->Write();
      h0J_qcd_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_anti->Write();
    }

    TDirectory *d0_qcd_taupt_iso =fout->mkdir("et_0jet_qcd_taupt_iso");
    d0_qcd_taupt_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_taupt_iso->SetName(name.c_str());
      h0LT_qcd_taupt_iso->Add(h0J_qcd_taupt_iso);
      h0LT_qcd_taupt_iso->Write();
    }
    else{
      h0LT_qcd_taupt_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_taupt_iso->Write();
      h0J_qcd_taupt_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_taupt_iso->Write();
    }
    TDirectory *d0_qcd_taupt_anti =fout->mkdir("et_0jet_qcd_taupt_anti");
    d0_qcd_taupt_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_taupt_anti->SetName(name.c_str());
      h0LT_qcd_taupt_anti->Add(h0J_qcd_taupt_anti);
      h0LT_qcd_taupt_anti->Write();
    }
    else{
      h0LT_qcd_taupt_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_taupt_anti->Write();
      h0J_qcd_taupt_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_taupt_anti->Write();
    }

    TDirectory *d0_qcd_pth_iso =fout->mkdir("et_0jet_qcd_pth_iso");
    d0_qcd_pth_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_pth_iso->SetName(name.c_str());
      h0LT_qcd_pth_iso->Add(h0J_qcd_pth_iso);
      h0LT_qcd_pth_iso->Write();
    }
    else{
      h0LT_qcd_pth_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_pth_iso->Write();
      h0J_qcd_pth_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_pth_iso->Write();
    }
    TDirectory *d0_qcd_pth_anti =fout->mkdir("et_0jet_qcd_pth_anti");
    d0_qcd_pth_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_pth_anti->SetName(name.c_str());
      h0LT_qcd_pth_anti->Add(h0J_qcd_pth_anti);
      h0LT_qcd_pth_anti->Write();
    }
    else{
      h0LT_qcd_pth_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_pth_anti->Write();
      h0J_qcd_pth_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_pth_anti->Write();
    }

    TDirectory *d0_qcd_mjj_iso =fout->mkdir("et_0jet_qcd_mjj_iso");
    d0_qcd_mjj_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_mjj_iso->SetName(name.c_str());
      h0LT_qcd_mjj_iso->Add(h0J_qcd_mjj_iso);
      h0LT_qcd_mjj_iso->Write();
    }
    else{
      h0LT_qcd_mjj_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_mjj_iso->Write();
      h0J_qcd_mjj_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_mjj_iso->Write();
    }
    TDirectory *d0_qcd_mjj_anti =fout->mkdir("et_0jet_qcd_mjj_anti");
    d0_qcd_mjj_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_mjj_anti->SetName(name.c_str());
      h0LT_qcd_mjj_anti->Add(h0J_qcd_mjj_anti);
      h0LT_qcd_mjj_anti->Write();
    }
    else{
      h0LT_qcd_mjj_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_mjj_anti->Write();
      h0J_qcd_mjj_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_mjj_anti->Write();
    }

    TDirectory *d0_qcd_cat0jetlow_iso =fout->mkdir("et_0jet_qcd_cat0jetlow_iso");
    d0_qcd_cat0jetlow_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_cat0jetlow_iso->SetName(name.c_str());
      h0LT_qcd_cat0jetlow_iso->Add(h0J_qcd_cat0jetlow_iso);
      h0LT_qcd_cat0jetlow_iso->Write();
    }
    else{
      h0LT_qcd_cat0jetlow_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_cat0jetlow_iso->Write();
      h0J_qcd_cat0jetlow_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_cat0jetlow_iso->Write();
    }
    TDirectory *d0_qcd_cat0jetlow_anti =fout->mkdir("et_0jet_qcd_cat0jetlow_anti");
    d0_qcd_cat0jetlow_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_cat0jetlow_anti->SetName(name.c_str());
      h0LT_qcd_cat0jetlow_anti->Add(h0J_qcd_cat0jetlow_anti);
      h0LT_qcd_cat0jetlow_anti->Write();
    }
    else{
      h0LT_qcd_cat0jetlow_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_cat0jetlow_anti->Write();
      h0J_qcd_cat0jetlow_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_cat0jetlow_anti->Write();
    }

    TDirectory *d0_qcd_cat0jethigh_iso =fout->mkdir("et_0jet_qcd_cat0jethigh_iso");
    d0_qcd_cat0jethigh_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_cat0jethigh_iso->SetName(name.c_str());
      h0LT_qcd_cat0jethigh_iso->Add(h0J_qcd_cat0jethigh_iso);
      h0LT_qcd_cat0jethigh_iso->Write();
    }
    else{
      h0LT_qcd_cat0jethigh_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_cat0jethigh_iso->Write();
      h0J_qcd_cat0jethigh_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_cat0jethigh_iso->Write();
    }
    TDirectory *d0_qcd_cat0jethigh_anti =fout->mkdir("et_0jet_qcd_cat0jethigh_anti");
    d0_qcd_cat0jethigh_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_cat0jethigh_anti->SetName(name.c_str());
      h0LT_qcd_cat0jethigh_anti->Add(h0J_qcd_cat0jethigh_anti);
      h0LT_qcd_cat0jethigh_anti->Write();
    }
    else{
      h0LT_qcd_cat0jethigh_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_cat0jethigh_anti->Write();
      h0J_qcd_cat0jethigh_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_cat0jethigh_anti->Write();
    }

    TDirectory *d0_qcd_catboosted1_iso =fout->mkdir("et_0jet_qcd_catboosted1_iso");
    d0_qcd_catboosted1_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_catboosted1_iso->SetName(name.c_str());
      h0LT_qcd_catboosted1_iso->Add(h0J_qcd_catboosted1_iso);
      h0LT_qcd_catboosted1_iso->Write();
    }
    else{
      h0LT_qcd_catboosted1_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_catboosted1_iso->Write();
      h0J_qcd_catboosted1_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_catboosted1_iso->Write();
    }
    TDirectory *d0_qcd_catboosted1_anti =fout->mkdir("et_0jet_qcd_catboosted1_anti");
    d0_qcd_catboosted1_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_catboosted1_anti->SetName(name.c_str());
      h0LT_qcd_catboosted1_anti->Add(h0J_qcd_catboosted1_anti);
      h0LT_qcd_catboosted1_anti->Write();
    }
    else{
      h0LT_qcd_catboosted1_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_catboosted1_anti->Write();
      h0J_qcd_catboosted1_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_catboosted1_anti->Write();
    }

    TDirectory *d0_qcd_catboosted2_iso =fout->mkdir("et_0jet_qcd_catboosted2_iso");
    d0_qcd_catboosted2_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_catboosted2_iso->SetName(name.c_str());
      h0LT_qcd_catboosted2_iso->Add(h0J_qcd_catboosted2_iso);
      h0LT_qcd_catboosted2_iso->Write();
    }
    else{
      h0LT_qcd_catboosted2_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_catboosted2_iso->Write();
      h0J_qcd_catboosted2_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_catboosted2_iso->Write();
    }
    TDirectory *d0_qcd_catboosted2_anti =fout->mkdir("et_0jet_qcd_catboosted2_anti");
    d0_qcd_catboosted2_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_catboosted2_anti->SetName(name.c_str());
      h0LT_qcd_catboosted2_anti->Add(h0J_qcd_catboosted2_anti);
      h0LT_qcd_catboosted2_anti->Write();
    }
    else{
      h0LT_qcd_catboosted2_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_catboosted2_anti->Write();
      h0J_qcd_catboosted2_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_catboosted2_anti->Write();
    }

    TDirectory *d0_qcd_catvbflow_iso =fout->mkdir("et_0jet_qcd_catvbflow_iso");
    d0_qcd_catvbflow_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_catvbflow_iso->SetName(name.c_str());
      h0LT_qcd_catvbflow_iso->Add(h0J_qcd_catvbflow_iso);
      h0LT_qcd_catvbflow_iso->Write();
    }
    else{
      h0LT_qcd_catvbflow_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_catvbflow_iso->Write();
      h0J_qcd_catvbflow_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_catvbflow_iso->Write();
    }
    TDirectory *d0_qcd_catvbflow_anti =fout->mkdir("et_0jet_qcd_catvbflow_anti");
    d0_qcd_catvbflow_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_catvbflow_anti->SetName(name.c_str());
      h0LT_qcd_catvbflow_anti->Add(h0J_qcd_catvbflow_anti);
      h0LT_qcd_catvbflow_anti->Write();
    }
    else{
      h0LT_qcd_catvbflow_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_catvbflow_anti->Write();
      h0J_qcd_catvbflow_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_catvbflow_anti->Write();
    }

    TDirectory *d0_qcd_catvbfhigh_iso =fout->mkdir("et_0jet_qcd_catvbfhigh_iso");
    d0_qcd_catvbfhigh_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_catvbfhigh_iso->SetName(name.c_str());
      h0LT_qcd_catvbfhigh_iso->Add(h0J_qcd_catvbfhigh_iso);
      h0LT_qcd_catvbfhigh_iso->Write();
    }
    else{
      h0LT_qcd_catvbfhigh_iso->SetName(name.c_str()+postfixLT);
      h0LT_qcd_catvbfhigh_iso->Write();
      h0J_qcd_catvbfhigh_iso->SetName(name.c_str()+postfixJ);
      h0J_qcd_catvbfhigh_iso->Write();
    }
    TDirectory *d0_qcd_catvbfhigh_anti =fout->mkdir("et_0jet_qcd_catvbfhigh_anti");
    d0_qcd_catvbfhigh_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_qcd_catvbfhigh_anti->SetName(name.c_str());
      h0LT_qcd_catvbfhigh_anti->Add(h0J_qcd_catvbfhigh_anti);
      h0LT_qcd_catvbfhigh_anti->Write();
    }
    else{
      h0LT_qcd_catvbfhigh_anti->SetName(name.c_str()+postfixLT);
      h0LT_qcd_catvbfhigh_anti->Write();
      h0J_qcd_catvbfhigh_anti->SetName(name.c_str()+postfixJ);
      h0J_qcd_catvbfhigh_anti->Write();
    }

    TDirectory *d0looseSS_qcd_iso =fout->mkdir("et_0SSloose_qcd_iso");
    d0looseSS_qcd_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0SSlooseLT_qcd_iso->SetName(name.c_str());
      h0SSlooseLT_qcd_iso->Add(h0SSlooseJ_qcd_iso);
      h0SSlooseLT_qcd_iso->Write();
    }
    else{
      h0SSlooseLT_qcd_iso->SetName(name.c_str()+postfixLT);
      h0SSlooseLT_qcd_iso->Write();
      h0SSlooseJ_qcd_iso->SetName(name.c_str()+postfixJ);
      h0SSlooseJ_qcd_iso->Write();
    }

    TDirectory *d0looseSS_qcd_anti =fout->mkdir("et_0SSloose_qcd_anti");
    d0looseSS_qcd_anti->cd();

    if (sample=="data_obs" or sample=="W"){
      h0SSlooseLT_qcd_anti->SetName(name.c_str());
      h0SSlooseLT_qcd_anti->Add(h0SSlooseJ_qcd_anti);
      h0SSlooseLT_qcd_anti->Write();
    }
    else{
      h0SSlooseLT_qcd_anti->SetName(name.c_str()+postfixLT);
      h0SSlooseLT_qcd_anti->Write();
      h0SSlooseJ_qcd_anti->SetName(name.c_str()+postfixJ);
      h0SSlooseJ_qcd_anti->Write();
    }

    TDirectory *d0_w_mvis_iso =fout->mkdir("et_0jet_w_mvis_iso");
    d0_w_mvis_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_mvis_iso->SetName(name.c_str());
      h0LT_w_mvis_iso->Add(h0J_w_mvis_iso);
      h0LT_w_mvis_iso->Write();
    }
    else{
      h0LT_w_mvis_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_mvis_iso->Write();
      h0J_w_mvis_iso->SetName(name.c_str()+postfixJ);
      h0J_w_mvis_iso->Write();
    }
    TDirectory *d0_w_mvis_anti =fout->mkdir("et_0jet_w_mvis_anti");
    d0_w_mvis_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_mvis_anti->SetName(name.c_str());
      h0LT_w_mvis_anti->Add(h0J_w_mvis_anti);
      h0LT_w_mvis_anti->Write();
    }
    else{
      h0LT_w_mvis_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_mvis_anti->Write();
      h0J_w_mvis_anti->SetName(name.c_str()+postfixJ);
      h0J_w_mvis_anti->Write();
    }

    TDirectory *d1_w_mvis_iso =fout->mkdir("et_1jet_w_mvis_iso");
    d1_w_mvis_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h1LT_w_mvis_iso->SetName(name.c_str());
      h1LT_w_mvis_iso->Add(h1J_w_mvis_iso);
      h1LT_w_mvis_iso->Write();
    }
    else{
      h1LT_w_mvis_iso->SetName(name.c_str()+postfixLT);
      h1LT_w_mvis_iso->Write();
      h1J_w_mvis_iso->SetName(name.c_str()+postfixJ);
      h1J_w_mvis_iso->Write();
    }
    TDirectory *d1_w_mvis_anti =fout->mkdir("et_1jet_w_mvis_anti");
    d1_w_mvis_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h1LT_w_mvis_anti->SetName(name.c_str());
      h1LT_w_mvis_anti->Add(h1J_w_mvis_anti);
      h1LT_w_mvis_anti->Write();
    }
    else{
      h1LT_w_mvis_anti->SetName(name.c_str()+postfixLT);
      h1LT_w_mvis_anti->Write();
      h1J_w_mvis_anti->SetName(name.c_str()+postfixJ);
      h1J_w_mvis_anti->Write();
    }

    TDirectory *d2_w_mvis_iso =fout->mkdir("et_2jet_w_mvis_iso");
    d2_w_mvis_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h2LT_w_mvis_iso->SetName(name.c_str());
      h2LT_w_mvis_iso->Add(h2J_w_mvis_iso);
      h2LT_w_mvis_iso->Write();
    }
    else{
      h2LT_w_mvis_iso->SetName(name.c_str()+postfixLT);
      h2LT_w_mvis_iso->Write();
      h2J_w_mvis_iso->SetName(name.c_str()+postfixJ);
      h2J_w_mvis_iso->Write();
    }
    TDirectory *d2_w_mvis_anti =fout->mkdir("et_2jet_w_mvis_anti");
    d2_w_mvis_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h2LT_w_mvis_anti->SetName(name.c_str());
      h2LT_w_mvis_anti->Add(h2J_w_mvis_anti);
      h2LT_w_mvis_anti->Write();
    }
    else{
      h2LT_w_mvis_anti->SetName(name.c_str()+postfixLT);
      h2LT_w_mvis_anti->Write();
      h2J_w_mvis_anti->SetName(name.c_str()+postfixJ);
      h2J_w_mvis_anti->Write();
    }


    TDirectory *d0_w_iso =fout->mkdir("et_0jet_w_iso");
    d0_w_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_iso->SetName(name.c_str());
      h0LT_w_iso->Add(h0J_w_iso);
      h0LT_w_iso->Write();
    }
    else{
      h0LT_w_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_iso->Write();
      h0J_w_iso->SetName(name.c_str()+postfixJ);
      h0J_w_iso->Write();
    }
    TDirectory *d0_w_anti =fout->mkdir("et_0jet_w_anti");
    d0_w_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_anti->SetName(name.c_str());
      h0LT_w_anti->Add(h0J_w_anti);
      h0LT_w_anti->Write();
    }
    else{
      h0LT_w_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_anti->Write();
      h0J_w_anti->SetName(name.c_str()+postfixJ);
      h0J_w_anti->Write();
    }




    TDirectory *d0_w_taupt_iso =fout->mkdir("et_0jet_w_taupt_iso");
    d0_w_taupt_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_taupt_iso->SetName(name.c_str());
      h0LT_w_taupt_iso->Add(h0J_w_taupt_iso);
      h0LT_w_taupt_iso->Write();
    }
    else{
      h0LT_w_taupt_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_taupt_iso->Write();
      h0J_w_taupt_iso->SetName(name.c_str()+postfixJ);
      h0J_w_taupt_iso->Write();
    }
    TDirectory *d0_w_taupt_anti =fout->mkdir("et_0jet_w_taupt_anti");
    d0_w_taupt_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_taupt_anti->SetName(name.c_str());
      h0LT_w_taupt_anti->Add(h0J_w_taupt_anti);
      h0LT_w_taupt_anti->Write();
    }
    else{
      h0LT_w_taupt_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_taupt_anti->Write();
      h0J_w_taupt_anti->SetName(name.c_str()+postfixJ);
      h0J_w_taupt_anti->Write();
    }

    TDirectory *d0_w_pth_iso =fout->mkdir("et_0jet_w_pth_iso");
    d0_w_pth_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_pth_iso->SetName(name.c_str());
      h0LT_w_pth_iso->Add(h0J_w_pth_iso);
      h0LT_w_pth_iso->Write();
    }
    else{
      h0LT_w_pth_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_pth_iso->Write();
      h0J_w_pth_iso->SetName(name.c_str()+postfixJ);
      h0J_w_pth_iso->Write();
    }
    TDirectory *d0_w_pth_anti =fout->mkdir("et_0jet_w_pth_anti");
    d0_w_pth_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_pth_anti->SetName(name.c_str());
      h0LT_w_pth_anti->Add(h0J_w_pth_anti);
      h0LT_w_pth_anti->Write();
    }
    else{
      h0LT_w_pth_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_pth_anti->Write();
      h0J_w_pth_anti->SetName(name.c_str()+postfixJ);
      h0J_w_pth_anti->Write();
    }

    TDirectory *d0_w_mjj_iso =fout->mkdir("et_0jet_w_mjj_iso");
    d0_w_mjj_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_mjj_iso->SetName(name.c_str());
      h0LT_w_mjj_iso->Add(h0J_w_mjj_iso);
      h0LT_w_mjj_iso->Write();
    }
    else{
      h0LT_w_mjj_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_mjj_iso->Write();
      h0J_w_mjj_iso->SetName(name.c_str()+postfixJ);
      h0J_w_mjj_iso->Write();
    }
    TDirectory *d0_w_mjj_anti =fout->mkdir("et_0jet_w_mjj_anti");
    d0_w_mjj_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_mjj_anti->SetName(name.c_str());
      h0LT_w_mjj_anti->Add(h0J_w_mjj_anti);
      h0LT_w_mjj_anti->Write();
    }
    else{
      h0LT_w_mjj_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_mjj_anti->Write();
      h0J_w_mjj_anti->SetName(name.c_str()+postfixJ);
      h0J_w_mjj_anti->Write();
    }

    TDirectory *d0_w_cat0jetlow_iso =fout->mkdir("et_0jet_w_cat0jetlow_iso");
    d0_w_cat0jetlow_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_cat0jetlow_iso->SetName(name.c_str());
      h0LT_w_cat0jetlow_iso->Add(h0J_w_cat0jetlow_iso);
      h0LT_w_cat0jetlow_iso->Write();
    }
    else{
      h0LT_w_cat0jetlow_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_cat0jetlow_iso->Write();
      h0J_w_cat0jetlow_iso->SetName(name.c_str()+postfixJ);
      h0J_w_cat0jetlow_iso->Write();
    }
    TDirectory *d0_w_cat0jetlow_anti =fout->mkdir("et_0jet_w_cat0jetlow_anti");
    d0_w_cat0jetlow_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_cat0jetlow_anti->SetName(name.c_str());
      h0LT_w_cat0jetlow_anti->Add(h0J_w_cat0jetlow_anti);
      h0LT_w_cat0jetlow_anti->Write();
    }
    else{
      h0LT_w_cat0jetlow_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_cat0jetlow_anti->Write();
      h0J_w_cat0jetlow_anti->SetName(name.c_str()+postfixJ);
      h0J_w_cat0jetlow_anti->Write();
    }

    TDirectory *d0_w_cat0jethigh_iso =fout->mkdir("et_0jet_w_cat0jethigh_iso");
    d0_w_cat0jethigh_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_cat0jethigh_iso->SetName(name.c_str());
      h0LT_w_cat0jethigh_iso->Add(h0J_w_cat0jethigh_iso);
      h0LT_w_cat0jethigh_iso->Write();
    }
    else{
      h0LT_w_cat0jethigh_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_cat0jethigh_iso->Write();
      h0J_w_cat0jethigh_iso->SetName(name.c_str()+postfixJ);
      h0J_w_cat0jethigh_iso->Write();
    }
    TDirectory *d0_w_cat0jethigh_anti =fout->mkdir("et_0jet_w_cat0jethigh_anti");
    d0_w_cat0jethigh_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_cat0jethigh_anti->SetName(name.c_str());
      h0LT_w_cat0jethigh_anti->Add(h0J_w_cat0jethigh_anti);
      h0LT_w_cat0jethigh_anti->Write();
    }
    else{
      h0LT_w_cat0jethigh_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_cat0jethigh_anti->Write();
      h0J_w_cat0jethigh_anti->SetName(name.c_str()+postfixJ);
      h0J_w_cat0jethigh_anti->Write();
    }

    TDirectory *d0_w_catboosted1_iso =fout->mkdir("et_0jet_w_catboosted1_iso");
    d0_w_catboosted1_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_catboosted1_iso->SetName(name.c_str());
      h0LT_w_catboosted1_iso->Add(h0J_w_catboosted1_iso);
      h0LT_w_catboosted1_iso->Write();
    }
    else{
      h0LT_w_catboosted1_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_catboosted1_iso->Write();
      h0J_w_catboosted1_iso->SetName(name.c_str()+postfixJ);
      h0J_w_catboosted1_iso->Write();
    }
    TDirectory *d0_w_catboosted1_anti =fout->mkdir("et_0jet_w_catboosted1_anti");
    d0_w_catboosted1_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_catboosted1_anti->SetName(name.c_str());
      h0LT_w_catboosted1_anti->Add(h0J_w_catboosted1_anti);
      h0LT_w_catboosted1_anti->Write();
    }
    else{
      h0LT_w_catboosted1_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_catboosted1_anti->Write();
      h0J_w_catboosted1_anti->SetName(name.c_str()+postfixJ);
      h0J_w_catboosted1_anti->Write();
    }

    TDirectory *d0_w_catboosted2_iso =fout->mkdir("et_0jet_w_catboosted2_iso");
    d0_w_catboosted2_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_catboosted2_iso->SetName(name.c_str());
      h0LT_w_catboosted2_iso->Add(h0J_w_catboosted2_iso);
      h0LT_w_catboosted2_iso->Write();
    }
    else{
      h0LT_w_catboosted2_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_catboosted2_iso->Write();
      h0J_w_catboosted2_iso->SetName(name.c_str()+postfixJ);
      h0J_w_catboosted2_iso->Write();
    }
    TDirectory *d0_w_catboosted2_anti =fout->mkdir("et_0jet_w_catboosted2_anti");
    d0_w_catboosted2_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_catboosted2_anti->SetName(name.c_str());
      h0LT_w_catboosted2_anti->Add(h0J_w_catboosted2_anti);
      h0LT_w_catboosted2_anti->Write();
    }
    else{
      h0LT_w_catboosted2_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_catboosted2_anti->Write();
      h0J_w_catboosted2_anti->SetName(name.c_str()+postfixJ);
      h0J_w_catboosted2_anti->Write();
    }

    TDirectory *d0_w_catvbflow_iso =fout->mkdir("et_0jet_w_catvbflow_iso");
    d0_w_catvbflow_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_catvbflow_iso->SetName(name.c_str());
      h0LT_w_catvbflow_iso->Add(h0J_w_catvbflow_iso);
      h0LT_w_catvbflow_iso->Write();
    }
    else{
      h0LT_w_catvbflow_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_catvbflow_iso->Write();
      h0J_w_catvbflow_iso->SetName(name.c_str()+postfixJ);
      h0J_w_catvbflow_iso->Write();
    }
    TDirectory *d0_w_catvbflow_anti =fout->mkdir("et_0jet_w_catvbflow_anti");
    d0_w_catvbflow_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_catvbflow_anti->SetName(name.c_str());
      h0LT_w_catvbflow_anti->Add(h0J_w_catvbflow_anti);
      h0LT_w_catvbflow_anti->Write();
    }
    else{
      h0LT_w_catvbflow_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_catvbflow_anti->Write();
      h0J_w_catvbflow_anti->SetName(name.c_str()+postfixJ);
      h0J_w_catvbflow_anti->Write();
    }

    TDirectory *d0_w_catvbfhigh_iso =fout->mkdir("et_0jet_w_catvbfhigh_iso");
    d0_w_catvbfhigh_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_catvbfhigh_iso->SetName(name.c_str());
      h0LT_w_catvbfhigh_iso->Add(h0J_w_catvbfhigh_iso);
      h0LT_w_catvbfhigh_iso->Write();
    }
    else{
      h0LT_w_catvbfhigh_iso->SetName(name.c_str()+postfixLT);
      h0LT_w_catvbfhigh_iso->Write();
      h0J_w_catvbfhigh_iso->SetName(name.c_str()+postfixJ);
      h0J_w_catvbfhigh_iso->Write();
    }
    TDirectory *d0_w_catvbfhigh_anti =fout->mkdir("et_0jet_w_catvbfhigh_anti");
    d0_w_catvbfhigh_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_w_catvbfhigh_anti->SetName(name.c_str());
      h0LT_w_catvbfhigh_anti->Add(h0J_w_catvbfhigh_anti);
      h0LT_w_catvbfhigh_anti->Write();
    }
    else{
      h0LT_w_catvbfhigh_anti->SetName(name.c_str()+postfixLT);
      h0LT_w_catvbfhigh_anti->Write();
      h0J_w_catvbfhigh_anti->SetName(name.c_str()+postfixJ);
      h0J_w_catvbfhigh_anti->Write();
    }


    TDirectory *dmt_w_iso =fout->mkdir("et_mt_w_iso");
    dmt_w_iso->cd();
    hmt_w_iso->SetName(name.c_str());
    hmt_w_iso->Write();

    TDirectory *dmt_w_anti =fout->mkdir("et_mt_w_anti");
    dmt_w_anti->cd();
    hmt_w_anti->SetName(name.c_str());
    hmt_w_anti->Write();

    TDirectory *d0_tt_iso =fout->mkdir("et_0jet_tt_iso");
    d0_tt_iso->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_tt_iso->SetName(name.c_str());
      h0LT_tt_iso->Add(h0J_tt_iso);
      h0LT_tt_iso->Write();
    }
    else{
      h0LT_tt_iso->SetName(name.c_str()+postfixLT);
      h0LT_tt_iso->Write();
      h0J_tt_iso->SetName(name.c_str()+postfixJ);
      h0J_tt_iso->Write();
    }

    TDirectory *d0_tt_anti =fout->mkdir("et_0jet_tt_anti");
    d0_tt_anti->cd();
    if (sample=="data_obs" or sample=="W"){
      h0LT_tt_anti->SetName(name.c_str());
      h0LT_tt_anti->Add(h0J_tt_anti);
      h0LT_tt_anti->Write();
    }
    else{
      h0LT_tt_anti->SetName(name.c_str()+postfixLT);
      h0LT_tt_anti->Write();
      h0J_tt_anti->SetName(name.c_str()+postfixJ);
      h0J_tt_anti->Write();
    }


    fout->Close();
    delete wmc;
} 


