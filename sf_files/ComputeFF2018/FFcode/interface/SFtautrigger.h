#include <TString.h> // Form

#include <iostream> // std::cerr, std::endl
#include <iomanip> 
#include <assert.h> // assert

const TH1* loadTH1(const TFile* inputFile, const std::string& histogramName)
{
  const TH1* histogram = dynamic_cast<TH1*>((const_cast<TFile*>(inputFile))->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = '" << histogramName << "' from input file !!" << std::endl;
    assert(0);
  }
  return histogram;
}

const TH2* loadTH2(const TFile* inputFile, const std::string& histogramName)
{
  const TH2* histogram = dynamic_cast<TH2*>((const_cast<TFile*>(inputFile))->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = '" << histogramName << "' from input file !!" << std::endl;
    assert(0);
  }
  return histogram;
}

double getEfficiency(double pt, double eta, double phi, const TH1* effHist, const TH2* etaPhi, const TH2* etaPhiAvg, int central_or_shift = 0)
{
  const TAxis* effHist_xAxis = effHist->GetXaxis();
  double ptMin = effHist_xAxis->GetXmin() + 1.e-1;
  double ptMax = effHist_xAxis->GetXmax() - 1.e-1;
  double pt_checked = pt;
  if ( pt_checked > ptMax ) pt_checked = ptMax;
  if ( pt_checked < ptMin ) pt_checked = ptMin;
  int effHist_idxBin = (const_cast<TH1*>(effHist))->FindBin(pt_checked);
  assert(effHist_idxBin >= 1 && effHist_idxBin <= effHist->GetNbinsX());
  double eff = effHist->GetBinContent(effHist_idxBin);

  double effErr = effHist->GetBinError(effHist_idxBin);
  if (central_or_shift==1) eff += effErr;
  if (central_or_shift==-1) eff -= effErr;

  // Adjust SF based on (eta, phi) location
  // keep eta barrel boundaries within SF region
  // but, for taus outside eta limits or with unralistic
  // phi values, return zero SF
  const TAxis* etaPhiAvg_xAxis = etaPhiAvg->GetXaxis();
  double etaMin = etaPhiAvg_xAxis->GetXmin() + 1.e-2;
  double etaMax = etaPhiAvg_xAxis->GetXmax() - 1.e-2;
  double eta_checked = eta;
  if ( eta_checked > etaMax ) eta_checked = etaMax;
  if ( eta_checked < etaMin ) eta_checked = etaMin;
  int etaPhiAvg_idxBinX = etaPhiAvg_xAxis->FindBin(eta_checked);
  assert(etaPhiAvg_idxBinX >= 1 && etaPhiAvg_idxBinX <= etaPhiAvg_xAxis->GetNbins());
  const TAxis* etaPhiAvg_yAxis = etaPhiAvg->GetYaxis();
  int etaPhiAvg_idxBinY = etaPhiAvg_yAxis->FindBin(phi);
  assert(etaPhiAvg_idxBinY >= 1 && etaPhiAvg_idxBinY <= etaPhiAvg_yAxis->GetNbins());
  double effCorr_etaPhi = etaPhi->GetBinContent((const_cast<TH2*>(etaPhi))->FindBin(eta_checked, phi));
  double effCorr_etaPhiAvg = etaPhiAvg->GetBinContent((const_cast<TH2*>(etaPhiAvg))->FindBin(eta_checked, phi));
  if ( effCorr_etaPhiAvg <= 0. ) {
    std::cerr << Form("One of the provided tau (eta, phi) values (%3.3f, %3.3f) is outside the boundary of triggering taus", eta, phi) << std::endl;
    std::cerr << "Returning efficiency = 0.0" << std::endl;
    return 0.;
  }  
  eff *= (effCorr_etaPhi/effCorr_etaPhiAvg);
  if ( eff > 1. ) eff = 1.;
  return eff;
}

double getDiTauEfficiencyData(double pt, double eta, double phi, int central_or_shift, TH1F* diTauData_, TH2F* diTauEtaPhiData_, TH2F* diTauEtaPhiAvgData_)
{
  return getEfficiency(pt, eta, phi, diTauData_, diTauEtaPhiData_, diTauEtaPhiAvgData_, central_or_shift);
}

double getDiTauEfficiencyMC(double pt, double eta, double phi, int central_or_shift, TH1F* diTauMC_, TH2F* diTauEtaPhiMC_, TH2F* diTauEtaPhiAvgMC_)
{
  return getEfficiency(pt, eta, phi, diTauMC_, diTauEtaPhiMC_, diTauEtaPhiAvgMC_, central_or_shift);
}

double getDiTauScaleFactor(double pt, double eta, double phi, int central_or_shift, TH1F* diTauMC_, TH2F* diTauEtaPhiMC_, TH2F* diTauEtaPhiAvgMC_,TH1F* diTauData_, TH2F* diTauEtaPhiData_, TH2F* diTauEtaPhiAvgData_)
{
  double effData = getDiTauEfficiencyData(pt, eta, phi, central_or_shift,diTauData_, diTauEtaPhiData_, diTauEtaPhiAvgData_);
  double effMC = getDiTauEfficiencyMC(pt, eta, phi, central_or_shift,diTauMC_, diTauEtaPhiMC_, diTauEtaPhiAvgMC_);
  if ( effMC < 1e-5 ) {
    std::cerr << "Eff MC is suspiciously low. Please contact Tau POG." << std::endl;
    std::cerr << Form(" - DiTau Trigger SF for Tau MVA:   pT: %3.3f   eta: %3.3f   phi: %3.3f", pt, eta, phi) << std::endl;
    std::cerr << Form(" - MC Efficiency = %3.3f", effMC) << std::endl;
    return 0.;
  }
  double sf = (effData/effMC);
  return sf;
}

