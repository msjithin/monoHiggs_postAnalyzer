void mutau_analyzer::fillHist( string histNumber , int muIndex, int tauIndex, bool isFakeBkg, float event_weight){
  string hNumber = histNumber;
  
  plotFill("muPt_"+hNumber,  my_muP4.Pt() , 30 , 20.0 , 80.0,  event_weight);
  plotFill("muEta_"+hNumber, my_muP4.Eta(), 48, -2.4, 2.4,  event_weight);
  plotFill("muPhi_"+hNumber, my_muP4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("muDz_"+hNumber,  muDz->at(muIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("muD0_"+hNumber,  muD0->at(muIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("muonID_"+hNumber, muIDbit->at(muIndex)>>1&1, 4, -2, 2,  event_weight); // electronID
  float relMuIso = ( muPFChIso->at(muIndex) + max( muPFNeuIso->at(muIndex) + muPFPhoIso->at(muIndex) - 0.5 *muPFPUIso->at(muIndex) , 0.0 )) / (my_muP4.Pt());
  plotFill("relMuIso_"+hNumber, relMuIso, 15, 0, 0.3,  event_weight);
  plotFill("muCharge_"+hNumber, muCharge->at(muIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  my_tauP4.Pt() , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, my_tauP4.Eta(), 25, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, my_tauP4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = my_muP4.DeltaR(my_tauP4);
  plotFill("deltaR_"+hNumber, deltaR , 40, 0, 6,  event_weight);
  double deltaPhi = DeltaPhi(muPhi->at(muIndex), tau_Phi->at(tauIndex));
  double deltaEta = fabs(muEta->at(muIndex) - tau_Eta->at(tauIndex));
  plotFill("deltaPhi_"+hNumber, deltaPhi , 30, -3.14, 3.14,  event_weight);
  plotFill("deltaEta_"+hNumber, deltaEta ,  25, -2.5, 2.5,  event_weight);
  
  plotFill("nVertex_"+hNumber, nVtx ,  24, 0, 60,  event_weight);
  plotFill("nJet_"+hNumber,  my_njets , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, my_metP4.Pt() , 20, 0, 200,  event_weight);
  plotFill("metLongXaxis_"+hNumber, my_metP4.Pt() , 10, 100, 200,  event_weight);
  plotFill("metPhi_"+hNumber, my_metP4.Phi() , 30, -3.14, 3.14,  event_weight);
  double mT_muMet = TMass_F( my_muP4.Pt(),my_muP4.Phi(), my_metP4.Pt(), my_metP4.Phi() );
  plotFill("mT_muMet_"+hNumber, mT_muMet , 20, 0, 200,event_weight);

  double visMass_mutau = (my_muP4+ my_tauP4).M();
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  double HiggsPt = (my_muP4+my_tauP4+my_metP4).Pt();
  plotFill("higgsPt_"+hNumber,HiggsPt ,  40, 0, 400,  event_weight);

  double tot_tr_mass = (my_muP4 + my_tauP4 + my_metP4 ).M();
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 16, 40, 200,  event_weight);
  if (tot_tr_mass >= 2000) tot_tr_mass = 1900;
  float TrMassBins[13]={ 40, 60, 90, 120, 150, 180, 210, 235, 260, 285, 325, 400, 2000};
  plotFill_customBinning("tot_TMass_full_"+hNumber, tot_tr_mass , 12, TrMassBins,  event_weight);

  int triggerBin1, triggerBin2, triggerBin3, triggerBin4;
  triggerBin1=triggerBin2=triggerBin3=triggerBin4=0;
  if( HLTEleMuX>>21&1 == 1 )  triggerBin3=3;
  if( HLTEleMuX>>60&1 == 1 ) triggerBin2=2;
  if( HLTTau>>0&1 == 1 )     triggerBin1=1;
  if(triggerBin1>0)
    plotFill("trigger_"+hNumber, triggerBin1 , 5, 0, 5,  event_weight);
  if(triggerBin2>0)
    plotFill("trigger_"+hNumber, triggerBin2 , 5, 0, 5,  event_weight);
  if(triggerBin3>0)
    plotFill("trigger_"+hNumber, triggerBin3 , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(is_MC){
    if(my_genmatching_l2==1) genMatchBin=1;
    else if(my_genmatching_l2==2) genMatchBin=2;
    else if(my_genmatching_l2==3) genMatchBin=3;
    else if(my_genmatching_l2==4) genMatchBin=4;
    else if(my_genmatching_l2==5) genMatchBin=5;
    else if(my_genmatching_l2==6) genMatchBin=6;
  }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
  
}

void mutau_analyzer::fillHist_nominal(string histNumber, float event_weight){
  string hNumber = histNumber;
  TLorentzVector tauP4, muP4, MET_P4;
  int muIndex, tauIndex;
  int genmatch1, genmatch2, njets;
  if (stage=="5_dyll")
    {
      muP4  = my_muP4_nom_5dyll;
      tauP4 = my_tauP4_nom_5dyll;
      MET_P4 = my_metP4_nom_5dyll;
      muIndex = MuIndex_nomdyll;
      tauIndex = TauIndex_nomdyll;
      genmatch1 = l1_genmatching_nomdyll;
      genmatch2 = l2_genmataching_nomdyll;
      njets= my_njets_nom_5dyll;
    }
  if (stage=="9_dyll")
    {
      muP4  = my_muP4_nom_9dyll;
      tauP4 = my_tauP4_nom_9dyll;
      MET_P4 = my_metP4_nom_9dyll;
      muIndex = MuIndex_nomdyll;
      tauIndex = TauIndex_nomdyll;
      genmatch1 = l1_genmatching_nomdyll;
      genmatch2 = l2_genmataching_nomdyll;
      njets= my_njets_nom_9dyll;
    }
  if (stage=="5")
    {
      muP4  = my_muP4_nom_5;
      tauP4 = my_tauP4_nom_5;
      MET_P4 = my_metP4_nom_5;
      muIndex = MuIndex_nom;
      tauIndex = TauIndex_nom;
      genmatch1 = l1_genmatching_nom;
      genmatch2 = l2_genmataching_nom;
      njets= my_njets_nom_5;
    }
  if (stage=="9")
    {
      muP4  = my_muP4_nom_9;
      tauP4 = my_tauP4_nom_9;
      MET_P4 = my_metP4_nom_9;
      muIndex = MuIndex_nom;
      tauIndex = TauIndex_nom;
      genmatch1 = l1_genmatching_nom;
      genmatch2 = l2_genmataching_nom;
      njets= my_njets_nom_9;
    }
  muIndex = MuIndex; tauIndex = TauIndex;
  plotFill("muPt_"+hNumber,  muP4.Pt() , 30 , 20.0 , 80.0,  event_weight);
  plotFill("muEta_"+hNumber, muP4.Eta(), 48, -2.4, 2.4,  event_weight);
  plotFill("muPhi_"+hNumber, muP4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("muDz_"+hNumber,  muDz->at(muIndex), 20, -0.2, 0.2,  event_weight);
  plotFill("muD0_"+hNumber,  muD0->at(muIndex), 48, -0.06, 0.06,  event_weight);
  plotFill("muonID_"+hNumber, muIDbit->at(muIndex)>>1&1, 4, -2, 2,  event_weight); // electronID
  float relMuIso = ( muPFChIso->at(muIndex) + max( muPFNeuIso->at(muIndex) + muPFPhoIso->at(muIndex) - 0.5 *muPFPUIso->at(muIndex) , 0.0 )) / (muP4.Pt());
  plotFill("relMuIso_"+hNumber, relMuIso, 15, 0, 0.3,  event_weight);
  plotFill("muCharge_"+hNumber, muCharge->at(muIndex), 8, -2, 2 ,  event_weight);
  
  plotFill("tauPt_"+hNumber,  tauP4.Pt() , 25 , 30 , 80,  event_weight);
  plotFill("tauEta_"+hNumber, tauP4.Eta(), 25, -2.5, 2.5,  event_weight);
  plotFill("tauPhi_"+hNumber, tauP4.Phi(), 30, -3.14, 3.14,  event_weight);
  plotFill("tauIso_"+hNumber, tau_byMediumDeepTau2017v2p1VSjet->at(tauIndex), 4, -2, 2,  event_weight);
  plotFill("tauDecayMode_"+hNumber, tau_DecayMode->at(tauIndex) , 12, 0, 12,  event_weight);
  plotFill("tauCharge_"+hNumber, tau_Charge->at(tauIndex), 8, -2, 2 ,  event_weight);
  plotFill("tauAntiEle_"+hNumber, tau_byTightDeepTau2017v2p1VSe->at(tauIndex), 8, -2, 2,  event_weight );
  plotFill("tauAntiMu_"+hNumber,  tau_byVLooseDeepTau2017v2p1VSmu->at(tauIndex), 8, -2, 2 ,  event_weight);
  double deltaR = muP4.DeltaR(tauP4);
  plotFill("deltaR_"+hNumber, deltaR , 40, 0, 6,  event_weight);
  double deltaPhi = DeltaPhi(muPhi->at(muIndex), tau_Phi->at(tauIndex));
  double deltaEta = fabs(muEta->at(muIndex) - tau_Eta->at(tauIndex));
  plotFill("deltaPhi_"+hNumber, deltaPhi , 30, -3.14, 3.14,  event_weight);
  plotFill("deltaEta_"+hNumber, deltaEta ,  25, -2.5, 2.5,  event_weight);
  
  plotFill("nVertex_"+hNumber, nVtx ,  24, 0, 60,  event_weight);
  plotFill("nJet_"+hNumber,  njets , 6, 0, 6,  event_weight);
  plotFill("met_"+hNumber, MET_P4.Pt() , 20, 0, 200,  event_weight);
  plotFill("metLongXaxis_"+hNumber, MET_P4.Pt() , 10, 100, 200,  event_weight);
  plotFill("metPhi_"+hNumber, MET_P4.Phi() , 30, -3.14, 3.14,  event_weight);
  double mT_muMet = TMass_F( muP4.Pt(),muP4.Phi(), MET_P4.Pt(), MET_P4.Phi() );
  plotFill("mT_muMet_"+hNumber, mT_muMet , 20, 0, 200,event_weight);

  double visMass_mutau = (muP4+ tauP4).M();
  plotFill("visMass_"+hNumber, visMass_mutau , 30, 50, 200,  event_weight);
  
  double HiggsPt = (muP4+tauP4+MET_P4).Pt();
  plotFill("higgsPt_"+hNumber,HiggsPt ,  40, 0, 400,  event_weight);

  double tot_tr_mass = (muP4 + tauP4 + MET_P4 ).M();
  plotFill("tot_TMass_"+hNumber, tot_tr_mass , 16, 40, 200,  event_weight);
  if (tot_tr_mass >= 2000) tot_tr_mass = 1900;
  float TrMassBins[13]={ 40, 60, 90, 120, 150, 180, 210, 235, 260, 285, 325, 400, 2000};
  plotFill_customBinning("tot_TMass_full_"+hNumber, tot_tr_mass , 12, TrMassBins,  event_weight);

  int triggerBin1, triggerBin2, triggerBin3, triggerBin4;
  triggerBin1=triggerBin2=triggerBin3=triggerBin4=0;
  if( HLTEleMuX>>21&1 == 1 )  triggerBin3=3;
  if( HLTEleMuX>>60&1 == 1 ) triggerBin2=2;
  if( HLTTau>>0&1 == 1 )     triggerBin1=1;
  if(triggerBin1>0)
    plotFill("trigger_"+hNumber, triggerBin1 , 5, 0, 5,  event_weight);
  if(triggerBin2>0)
    plotFill("trigger_"+hNumber, triggerBin2 , 5, 0, 5,  event_weight);
  if(triggerBin3>0)
    plotFill("trigger_"+hNumber, triggerBin3 , 5, 0, 5,  event_weight);
  
  int genMatchBin=0;
  if(is_MC){
    if(genmatch2==1) genMatchBin=1;
    else if(genmatch2==2) genMatchBin=2;
    else if(genmatch2==3) genMatchBin=3;
    else if(genmatch2==4) genMatchBin=4;
    else if(genmatch2==5) genMatchBin=5;
    else if(genmatch2==6) genMatchBin=6;
  }
  plotFill("genMatch_"+hNumber, genMatchBin ,7, 0, 7,  event_weight);
  //if(debug)cout <<"plots filled for "<<hNumber<<endl;
}
