void etau_analyzer::fillHist( string histNumber , int eleIndex, int tauIndex, float event_weight, bool isMC){
  string hNumber = histNumber;
  plotFill("elePt_"+hNumber,  elePt->at(eleIndex) , 38 , 24 , 100,  event_weight);
  plotFill("eleEta_"+hNumber, eleEta->at(eleIndex), 48, -2.4, 2.4,  event_weight);
  plotFill("elePhi_"+hNumber, elePhi->at(eleIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("eleDz_"+hNumber,  eleDz->at(eleIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("eleD0_"+hNumber,  eleD0->at(eleIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("electronID_"+hNumber, eleIDbit->at(eleIndex)>>8&1, 4, -2, 2,  event_weight); // electronID
  float relEleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (elePt->at(eleIndex));
  plotFill("relEleIso_"+hNumber, relEleIso, 15, 0, 0.3,  event_weight);
  plotFill("eleCharge_"+hNumber, eleCharge->at(eleIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 45, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), tau_Phi->at(tauIndex), tau_Eta->at(tauIndex));
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(eleIndex, tauIndex, -1);
  plotFill("nJet_"+hNumber,  jetCand.size() , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, pfMET , 40, 0, 200,  event_weight);
  
  double mT_eleMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
  plotFill("mT_eleMet_"+hNumber, mT_eleMet , 30, 0, 150,  event_weight);

  TLorentzVector myTau; 
  myTau.SetPtEtaPhiE(tau_Pt->at(tauIndex),tau_Eta->at(tauIndex),tau_Phi->at(tauIndex), tau_Energy->at(tauIndex));
  TLorentzVector myEle; 
  myEle.SetPtEtaPhiE(elePt->at(eleIndex), eleEta->at(eleIndex), elePhi->at(eleIndex), eleE->at(eleIndex));
  double visMass_mutau = VisMass_F(myTau, myEle);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  double HiggsPt = pTvecsum_F(myEle, myTau, myMet );
  plotFill("higgsPt_"+hNumber,HiggsPt , 25, 0, 250,  event_weight);

  double tot_tr_mass = TotTMass_F(myEle, myTau, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  int triggerBin=0;
  if( HLTTau>>1&1 == 1 )     triggerBin=4;
  else if( HLTEleMuX>>3&1 == 1 ) triggerBin=1;
  else if( HLTEleMuX>>61&1 == 1 ) triggerBin=2;
  else if( HLTEleMuX>>5&1 == 1 ) triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(isMC){
    if(found_GenMatch(1)) genMatchBin=1;
    else if(found_GenMatch(2)) genMatchBin=2;
    else if(found_GenMatch(3)) genMatchBin=3;
    else if(found_GenMatch(4)) genMatchBin=4;
    else if(found_GenMatch(5)) genMatchBin=5;
    else if(found_GenMatch(6)) genMatchBin=6;
  }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}
void etau_analyzer::fillHist( string histNumber , TLorentzVector eleP4, TLorentzVector tauP4, int eleIndex, int tauIndex, float event_weight, bool isMC){
  string hNumber = histNumber;
  plotFill("elePt_"+hNumber,  elePt->at(eleIndex) , 38 , 24 , 100,  event_weight);
  plotFill("eleEta_"+hNumber, eleEta->at(eleIndex), 48, -2.4, 2.4,  event_weight);
  plotFill("elePhi_"+hNumber, elePhi->at(eleIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("eleDz_"+hNumber,  eleDz->at(eleIndex), 40, -0.5, 0.5,  event_weight);
  plotFill("eleD0_"+hNumber,  eleD0->at(eleIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("electronID_"+hNumber, eleIDbit->at(eleIndex)>>8&1, 4, -2, 2,  event_weight); 
  float relEleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (elePt->at(eleIndex));
  plotFill("relEleIso_"+hNumber, relEleIso, 15, 0, 0.3,  event_weight);
  plotFill("eleCharge_"+hNumber, eleCharge->at(eleIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 45, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), tau_Phi->at(tauIndex), tau_Eta->at(tauIndex));
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);

  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(eleIndex, tauIndex, -1);
  plotFill("nJet_"+hNumber, jetCand.size() , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, pfMET , 40, 0, 200,  event_weight);

  double mT_muMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
  plotFill("mT_eleMet_"+hNumber, mT_muMet , 30, 0, 150,  event_weight);
  
  double visMass_mutau = VisMass_F(eleP4, tauP4);
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  double HiggsPt = pTvecsum_F(eleP4, tauP4, myMet);
  plotFill("higgsPt_"+hNumber,HiggsPt , 25, 0, 250,  event_weight);

  double tot_tr_mass = TotTMass_F(eleP4, tauP4, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  int triggerBin=0;
  if( HLTTau>>1&1 == 1 )     triggerBin=4;
  else if( HLTEleMuX>>3&1 == 1 ) triggerBin=1;
  else if( HLTEleMuX>>61&1 == 1 ) triggerBin=2;
  else if( HLTEleMuX>>5&1 == 1 ) triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(isMC){
    if(found_GenMatch(1)) genMatchBin=1;
    else if(found_GenMatch(2)) genMatchBin=2;
    else if(found_GenMatch(3)) genMatchBin=3;
    else if(found_GenMatch(4)) genMatchBin=4;
    else if(found_GenMatch(5)) genMatchBin=5;
    else if(found_GenMatch(6)) genMatchBin=6;
  }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
}
void etau_analyzer::fillHist_dyll( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight, bool isMC){
  string hNumber = histNumber;
  plotFill("elePt_"+hNumber,  elePt->at(eleIndex) , 38 , 24 , 100,  event_weight);
  plotFill("eleEta_"+hNumber, eleEta->at(eleIndex), 48, -2.4, 2.4,  event_weight);
  plotFill("elePhi_"+hNumber, elePhi->at(eleIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("eleDz_"+hNumber,  eleDz->at(eleIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("eleD0_"+hNumber,  eleD0->at(eleIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("electronID_"+hNumber, eleIDbit->at(eleIndex)>>8&1, 4, -2, 2,  event_weight); // electronID
  float relEleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (elePt->at(eleIndex));
  plotFill("relEleIso_"+hNumber, relEleIso, 15, 0, 0.3,  event_weight);
  plotFill("eleCharge_"+hNumber, eleCharge->at(eleIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tau_Pt->at(tauIndex) , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tau_Eta->at(tauIndex), 45, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tau_Phi->at(tauIndex), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  
  double deltaR = delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), elePhi->at(ele2Index), eleEta->at(ele2Index));
  plotFill("deltaR_"+hNumber, deltaR , 30, 0, 6,  event_weight);
  std::vector<int> jetCand;       jetCand.clear();
  jetCand=getJetCand(eleIndex, -1, ele2Index);
  plotFill("nJet_"+hNumber,  jetCand.size() , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, pfMET , 40, 0, 200,  event_weight);
  
  double mT_eleMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
  plotFill("mT_eleMet_"+hNumber, mT_eleMet , 30, 0, 150,  event_weight);

  TLorentzVector myTau; 
  myTau.SetPtEtaPhiE(tau_Pt->at(tauIndex),tau_Eta->at(tauIndex),tau_Phi->at(tauIndex), tau_Energy->at(tauIndex));
  TLorentzVector myEle; 
  myEle.SetPtEtaPhiE(elePt->at(eleIndex), eleEta->at(eleIndex), elePhi->at(eleIndex), eleE->at(eleIndex));
  TLorentzVector myEle2;
  myEle2.SetPtEtaPhiE(elePt->at(ele2Index), eleEta->at(ele2Index), elePhi->at(ele2Index), eleE->at(ele2Index));
  double visMass_eletau = VisMass_F(myEle2, myEle);
  plotFill("visMass_"+hNumber, visMass_eletau , 30, 50, 200,  event_weight);
  
  TLorentzVector myMet;
  myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
  double HiggsPt = pTvecsum_F(myEle, myEle2, myMet );
  plotFill("higgsPt_"+hNumber,HiggsPt , 25, 0, 250,  event_weight);

  double tot_tr_mass = TotTMass_F(myEle, myEle2, myMet );
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 20, 0, 200,  event_weight);

  int triggerBin=0;
  if( HLTTau>>1&1 == 1 )     triggerBin=4;
  else if( HLTEleMuX>>3&1 == 1 ) triggerBin=1;
  else if( HLTEleMuX>>61&1 == 1 ) triggerBin=2;
  else if( HLTEleMuX>>5&1 == 1 ) triggerBin=3;
  plotFill("trigger_"+hNumber, triggerBin , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(isMC)
    {
      if(found_GenMatch(1)) genMatchBin=1;
      else if(found_GenMatch(2)) genMatchBin=2;
      else if(found_GenMatch(3)) genMatchBin=3;
      else if(found_GenMatch(4)) genMatchBin=4;
      else if(found_GenMatch(5)) genMatchBin=5;
      else if(found_GenMatch(6)) genMatchBin=6;
    }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);  
}

void etau_analyzer::makeTestPlot( string histNumber , int eleIndex, int ele2Index, int tauIndex, float event_weight){
  string hNumber = histNumber;
  std::vector<int> tmpCand; tmpCand.clear();
  for(int iEle=0;iEle<nEle;iEle++)
    {
      tmpCand.push_back(iEle);
    }
  plotFill("elePt_"+hNumber,  elePt->at(tmpCand[0]) , 38 , 24 , 100,  event_weight);
  //cout<<"     elePt_"<<hNumber<<" = "<< elePt->at(tmpCand[0])<<endl;
}
