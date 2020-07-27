double etau_analyzer::getScaleFactors( int eleIndex , int tauIndex, bool fakeBkg , bool isMC, bool debug)
{
  double rv_sf=1.0;
  double eleRecoSF_corr=1.0;
  double eleEffSF_corr=1.0;
  double eletrgsf_tmp=1.0;
  double eletrgsf=1.0;
  double eletrgsf_ele32=1.0;
  double eletrgsf_ele35=1.0;
  double eletrgsf_ele27=1.0;
  double sf_Zvtx=0.991;
  double sf_tauidSF_m = 1.0;
  double sf_tauidSF_vvvl = 1.0;
  double sf_tauesSF = 1.0;
  double sf_fakeEle = 1.0; double sf_fakeMu = 1.0;
  double sf_taufesSF = 1.0;
  double sf_tauTrg = 1.0;
  //int genMatchTau = 0;
  //if(isMC) genMatchTau = gen_matching();
  
  eleRecoSF_corr=h_eleRecoSF_highpt->GetBinContent(h_eleRecoSF_highpt->GetXaxis()->FindBin(eleSCEta->at(eleIndex)),h_eleRecoSF_highpt->GetYaxis()->FindBin(elePt->at(eleIndex)));
  if (debug==true ) std::cout<<"eleRecoSF_corr =  "<< eleRecoSF_corr<<std::endl;
  eleEffSF_corr=h_eleIDSF->GetBinContent(h_eleIDSF->GetXaxis()->FindBin(eleSCEta->at(eleIndex)),h_eleIDSF->GetYaxis()->FindBin(elePt->at(eleIndex)));
  if (debug==true ) std::cout<<"eleEffSF_corr =  "<< eleEffSF_corr<<std::endl;
  if (debug==true ) std::cout<<"This works line 269 "<<std::endl;
  //eletrgsf = EletriggerSF(elePt->at(eleIndex), eleEta->at(eleIndex));
  eletrgsf_tmp =h_eleTrgSF_1->GetBinContent(h_eleTrgSF_1->GetXaxis()->FindBin(eleSCEta->at(eleIndex)),h_eleTrgSF_1->GetYaxis()->FindBin(elePt->at(eleIndex)));
  if(eletrgsf_tmp>0)
    eletrgsf=eletrgsf_tmp;
  //eletrgsf =h_eleTrgSF_2->GetBinContent(h_eleTrgSF_2->GetXaxis()->FindBin(eleSCEta->at(eleIndex)),h_eleTrgSF_2->GetYaxis()->FindBin(elePt->at(eleIndex)));
  if (debug==true ) std::cout<<"eletrgsf =  "<< eletrgsf << "    Line 1168"<<std::endl;
  


  if( found_GenMatch(5))
    {
      sf_tauidSF_m = h_tauidSF_m->GetBinContent(h_tauidSF_m->GetXaxis()->FindBin(tau_DecayMode->at(tauIndex)));
      //sf_tauidSF_m = fn_tauIDSF_m->Eval(tau_Pt->at(reco_tau[0]));
      sf_tauidSF_vvvl = h_tauidSF_vvvl->GetBinContent(h_tauidSF_vvvl->GetXaxis()->FindBin(tau_DecayMode->at(tauIndex)));
      //sf_tauidSF_vvvl = fn_tauIDSF_vvl->Eval(tau_Pt->at(reco_tau[0]));
      sf_tauesSF = h_tauesSF->GetBinContent(h_tauesSF->GetXaxis()->FindBin(tau_DecayMode->at(tauIndex)));
    }
  //if(genMatchTau==1 || genMatchTau==3)
  //sf_fakeEle = h_tauFakeEleSF->GetBinContent(h_tauFakeEleSF->GetXaxis()->FindBin(abs(tau_Eta->at(tauIndex))));
  //if(genMatchTau==2 || genMatchTau==4)
  //sf_fakeMu = h_tauFakeMuSF->GetBinContent(h_tauFakeMuSF->GetXaxis()->FindBin(abs(tau_Eta->at(tauIndex))));
  if(found_GenMatch(2) || found_GenMatch(4)){
    if(tau_DecayMode->at(tauIndex)==0)
      {
	if(abs(tau_Eta->at(tauIndex)) < 0.4 ) sf_fakeMu=1.14;
	if(abs(tau_Eta->at(tauIndex)) > 0.4 
	   && abs(tau_Eta->at(tauIndex)) < 0.8 ) sf_fakeMu=1.0;
	if(abs(tau_Eta->at(tauIndex)) > 0.8
           && abs(tau_Eta->at(tauIndex)) < 1.2 ) sf_fakeMu=0.87;
	if(abs(tau_Eta->at(tauIndex)) > 1.2
           && abs(tau_Eta->at(tauIndex)) < 1.7 ) sf_fakeMu=0.52;
	if(abs(tau_Eta->at(tauIndex)) > 1.7
           && abs(tau_Eta->at(tauIndex)) < 2.3 ) sf_fakeMu=1.47;
      }
    if(tau_DecayMode->at(tauIndex)==1)
      {
	if(abs(tau_Eta->at(tauIndex)) > 0.0
           && abs(tau_Eta->at(tauIndex)) < 0.4 ) sf_fakeMu=0.69;
      }
  }
  if(found_GenMatch(1) || found_GenMatch(3)){
    if(tau_DecayMode->at(tauIndex)==0)
      {
	if(abs(tau_Eta->at(tauIndex)) < 1.479 ) sf_fakeEle=0.98;
	if(abs(tau_Eta->at(tauIndex)) > 1.479 ) sf_fakeEle=0.80;
      }
    if(tau_DecayMode->at(tauIndex)==1)
      {
	if(abs(tau_Eta->at(tauIndex)) < 1.479 ) sf_fakeEle=1.07;
        if(abs(tau_Eta->at(tauIndex)) > 1.479 ) sf_fakeEle=0.64;
      }
  }
  if(tau_DecayMode->at(tauIndex)==0 && abs(tau_Eta->at(tauIndex))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(1);
  if(tau_DecayMode->at(tauIndex)==0 && abs(tau_Eta->at(tauIndex))>1.4 )  sf_taufesSF = h_taufesSF->Eval(3);
  if(tau_DecayMode->at(tauIndex)==1 && abs(tau_Eta->at(tauIndex))<=1.4 ) sf_taufesSF = h_taufesSF->Eval(5);
  if(tau_DecayMode->at(tauIndex)==1 && abs(tau_Eta->at(tauIndex))>1.4 )  sf_taufesSF = h_taufesSF->Eval(7);
  
  /* cout<<"tau trg sf = "<<tauSFs.getTriggerScaleFactor(  tau_Pt->at(tauIndex) */
  /* 						      , tau_Eta->at(tauIndex) */
  /* 						      , tau_Phi->at(tauIndex) */
  /* 						      , tau_DecayMode->at(tauIndex) */
  /* 						      )<<endl; */
  double tauPtCheck=0.0;
  if(tau_Pt->at(tauIndex) > 450 ) tauPtCheck = 450;
  else if ( tau_Pt->at(tauIndex) < 20 )  tauPtCheck = 20;
  
  if(tau_DecayMode->at(tauIndex)==0)  sf_tauTrg= h_tauTrgSF_dm0->GetBinContent(h_tauTrgSF_dm0->GetXaxis()->FindBin(tauPtCheck));
  if(tau_DecayMode->at(tauIndex)==1)  sf_tauTrg= h_tauTrgSF_dm1->GetBinContent(h_tauTrgSF_dm1->GetXaxis()->FindBin(tauPtCheck));
  if(tau_DecayMode->at(tauIndex)==10) sf_tauTrg= h_tauTrgSF_dm10->GetBinContent(h_tauTrgSF_dm10->GetXaxis()->FindBin(tauPtCheck));
  if(tau_DecayMode->at(tauIndex)==11) sf_tauTrg= h_tauTrgSF_dm11->GetBinContent(h_tauTrgSF_dm11->GetXaxis()->FindBin(tauPtCheck));
  //cout<<"sf_tauTrg = "<<sf_tauTrg<<endl;
  //event_weight=event_weight * sf_muID * sf_IsoEff * sf_muTrg * sf_tauidSF_m * sf_tauesSF * sf_fakeEle * (sf_fakeMu);
  if(fakeBkg)
    rv_sf = eleRecoSF_corr * eleEffSF_corr * eletrgsf * sf_Zvtx * sf_tauidSF_m * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF * sf_tauTrg;
  else
    rv_sf = eleRecoSF_corr * eleEffSF_corr * eletrgsf * sf_Zvtx * sf_tauidSF_m * sf_tauidSF_vvvl * sf_tauesSF * sf_fakeEle * sf_fakeMu * sf_taufesSF * sf_tauTrg;
  
  return rv_sf;
  //return 1.0;
}


float etau_analyzer::EletriggerSF(float pt, float eta){
   float sf = 1.0;
   if(fabs(eta) >= 0.0   && fabs(eta) < 0.8)
     {
       if(pt > 40.0  && pt < 50) sf = 0.79;
       if(pt > 50.0  && pt < 60) sf = 0.82;
       if(pt > 60.0  && pt < 100) sf = 0.85;
       if(pt > 100.0  && pt < 150) sf = 0.87;
       if(pt > 150.0  && pt < 200) sf = 0.88;
       if(pt > 200) sf = 0.89;
     }
   if(fabs(eta) >= 0.8   && fabs(eta) < 1.442 ) 
     {
       if(pt > 40.0  && pt < 50) sf = 0.77;
       if(pt > 50.0  && pt < 60) sf = 0.81;
       if(pt > 60.0  && pt < 100) sf = 0.85;
       if(pt > 100.0  && pt < 150) sf = 0.87;
       if(pt > 150.0  && pt < 300) sf = 0.89;
       if(pt > 300.0) sf = 0.87;
     }
   if(fabs(eta) >= 1.442   && fabs(eta) < 1.56 ) 
     {
       if(pt > 40.0  && pt < 50) sf = 0.73;
       if(pt > 50.0  && pt < 60) sf = 0.75;
       if(pt > 60.0  && pt < 100) sf = 0.76;
       if(pt > 100.0  && pt < 150) sf = 0.72;
       if(pt > 150.0  && pt < 300) sf = 0.78;
       if(pt > 300.0) sf = 0.67;
     }
   if(fabs(eta) >= 1.56   && fabs(eta) < 2.0 ) 
     {
       if(pt > 40.0  && pt < 50) sf = 0.80;
       if(pt > 50.0  && pt < 60) sf = 0.84;
       if(pt > 60.0  && pt < 100) sf = 0.87;
       if(pt > 100.0  && pt < 150) sf = 0.88;
       if(pt > 150.0  && pt < 300) sf = 0.89;
       if(pt > 300.0) sf = 0.87;
     }
   if(fabs(eta) >= 2.0   && fabs(eta) < 2.5 ) 
     {
       if(pt > 40.0  && pt < 50) sf = 0.73;
       if(pt > 50.0  && pt < 60) sf = 0.78;
       if(pt > 60.0  && pt < 100) sf = 0.83;
       if(pt > 100.0  && pt < 150) sf = 0.86;
       if(pt > 150.0  && pt < 300) sf = 0.89;
       if(pt > 300.0) sf = 0.86;
     }
   return sf;
}
