TLorentzVector etau_analyzer::metSysUnc(string uncType, TLorentzVector eventMetP4){
  
  TLorentzVector mymetvector = eventMetP4;
  //TLorentzVector mymetvector;
  //mymetvector.SetPtEtaPhiE(pfMET,0 ,pfMETPhi ,pfMET);
  TLorentzVector BosonP4, nuP4, nuP4tmp;
  TLorentzVector nu1P4, gentau1P4;
  TLorentzVector nu2P4, gentau2P4;
  TLorentzVector visGenP4;
  for(int i=0; i<nMC; i++)
    {
      if(mcPID->at(i)==23)
	BosonP4.SetPtEtaPhiE(mcPt->at(i), mcEta->at(i) , mcPhi->at(i) , mcE->at(i) );
    }
  //visGenP4=BosonP4;
  if(BosonP4.Pt()==0)
    {
      for(int i=0; i<nMC; i++)
	{
	  if(mcPID->at(i)==15)
	    gentau1P4.SetPtEtaPhiE(mcPt->at(i), mcEta->at(i) , mcPhi->at(i) , mcE->at(i) );
	  if(mcPID->at(i)==-15)
	    gentau2P4.SetPtEtaPhiE(mcPt->at(i), mcEta->at(i) , mcPhi->at(i) , mcE->at(i) );
	}
      BosonP4=gentau1P4+gentau2P4;
    }
  visGenP4=BosonP4;
  for(int i=0; i<nMC; i++)
    {
      if(abs(mcPID->at(i))==16 || abs(mcPID->at(i))==14 || abs(mcPID->at(i))==12)
	{
	  nuP4tmp.SetPtEtaPhiE(mcPt->at(i), mcEta->at(i) , mcPhi->at(i) , mcE->at(i) );
	  visGenP4=visGenP4-nuP4tmp;
	}
    }
  // apply recoil corrections on event-by-event basis (Type I PF MET)
  float _pfmet=mymetvector.Pt(); float _pfmetPhi=mymetvector.Phi();
  float pfmetcorr_ex_responseUp=_pfmet*cos(_pfmetPhi);     float pfmetcorr_ey_responseUp=_pfmet*sin(_pfmetPhi);
  float pfmetcorr_ex_responseDown=_pfmet*cos(_pfmetPhi);   float pfmetcorr_ey_responseDown=_pfmet*sin(_pfmetPhi);
  float pfmetcorr_ex_resolutionUp=_pfmet*cos(_pfmetPhi);   float pfmetcorr_ey_resolutionUp=_pfmet*sin(_pfmetPhi);
  float pfmetcorr_ex_resolutionDown=_pfmet*cos(_pfmetPhi); float pfmetcorr_ey_resolutionDown=_pfmet*sin(_pfmetPhi);
  float genpX = BosonP4.Px();
  float genpY = BosonP4.Py();
  float vispX = visGenP4.Px();
  float vispY = visGenP4.Py();
  int recoiljets = my_njets;
  int Process= MEtSys::ProcessType::BOSON;
  metSys.ApplyMEtSys(mymetvector.Px(),mymetvector.Py(),
		     genpX,genpY,
		     vispX,vispY,
		     recoiljets,
		     Process,
		     MEtSys::SysType::Response,
		     MEtSys::SysShift::Up,
		     pfmetcorr_ex_responseUp,pfmetcorr_ey_responseUp
		     );

  metSys.ApplyMEtSys(mymetvector.Px(),mymetvector.Py(),
		     genpX,genpY,
		     vispX,vispY,
		     recoiljets,
		     Process,
		     MEtSys::SysType::Response,
		     MEtSys::SysShift::Down,
		     pfmetcorr_ex_responseDown,pfmetcorr_ey_responseDown
		     );

  metSys.ApplyMEtSys(mymetvector.Px(),mymetvector.Py(),
		     genpX,genpY,
		     vispX,vispY,
		     recoiljets,
		     Process,
		     MEtSys::SysType::Resolution,
		     MEtSys::SysShift::Up,
		     pfmetcorr_ex_resolutionUp,pfmetcorr_ey_resolutionUp
		     );
  metSys.ApplyMEtSys(mymetvector.Px(),mymetvector.Py(),
		     genpX,genpY,
		     vispX,vispY,
		     recoiljets,
		     Process,
		     MEtSys::SysType::Resolution,
		     MEtSys::SysShift::Down,
		     pfmetcorr_ex_resolutionDown,pfmetcorr_ey_resolutionDown
		     );
  // cout<<" ********************************************"<<endl;
  // cout<<"nominal    "<<" "<<mymetvector.Px()<<"  "<<mymetvector.Py()<<endl;
  // cout <<"response"<<"\n"<<" "<<pfmetcorr_ex_responseUp<<" "<<pfmetcorr_ex_responseDown<<endl;
  // cout <<"response"<<"\n"<<" "<<pfmetcorr_ey_responseUp<<" "<<pfmetcorr_ey_responseDown<<endl;
  // cout <<"resolution"<<"\n"<<" "<<pfmetcorr_ex_resolutionUp<<" "<<pfmetcorr_ex_resolutionDown<<endl;
  // cout <<"resolution"<<"\n"<<" "<<pfmetcorr_ey_resolutionUp<<" "<<pfmetcorr_ey_resolutionDown<<endl;
  string njetName ="";
  // if (recoiljets ==0) njetName="0";
  // else if (recoiljets ==1) njetName="1";
  // else if (recoiljets >=2) njetName="2";
  TLorentzVector mymet_responseUp, mymet_responseDown, mymet_resolutionUp, mymet_resolutionDown ;
  if (unc_shift == "up" && uncType=="response")
    {
      mymet_responseUp.SetPxPyPzE(pfmetcorr_ex_responseUp,pfmetcorr_ey_responseUp,0,sqrt(pfmetcorr_ex_responseUp*pfmetcorr_ex_responseUp+pfmetcorr_ey_responseUp*pfmetcorr_ey_responseUp));
      //if(make_met_plot==true) { cout<<"filling met_response_up"<<endl; plotFill("met_response_up", mymet_responseUp.Pt(), 20, 0, 200,  1.0);}
      return mymet_responseUp;
    }
  else if (unc_shift == "up" && uncType=="resolution")
    {
      mymet_resolutionUp.SetPxPyPzE(pfmetcorr_ex_resolutionUp,pfmetcorr_ey_resolutionUp,0,sqrt(pfmetcorr_ex_resolutionUp*pfmetcorr_ex_resolutionUp+pfmetcorr_ey_resolutionUp*pfmetcorr_ey_resolutionUp));
      //if(make_met_plot==true) { cout<<"filling met_resolution_up"<<endl; plotFill("met_resolution_up", mymet_resolutionUp.Pt(), 20, 0, 200,  1.0);}
      return mymet_resolutionUp;
    }
  else if (unc_shift == "down" && uncType=="response")
    {
      mymet_responseDown.SetPxPyPzE(pfmetcorr_ex_responseDown,pfmetcorr_ey_responseDown,0,sqrt(pfmetcorr_ex_responseDown*pfmetcorr_ex_responseDown+pfmetcorr_ey_responseDown*pfmetcorr_ey_responseDown));
      //if(make_met_plot==true) { cout<<"filling met_response_down"<<endl; plotFill("met_response_down", mymet_responseDown.Pt(), 20, 0, 200,  1.0);}
      return mymet_responseDown;
    }
  else if (unc_shift == "down" && uncType=="resolution")
    {
      mymet_resolutionDown.SetPxPyPzE(pfmetcorr_ex_resolutionDown,pfmetcorr_ey_resolutionDown,0,sqrt(pfmetcorr_ex_resolutionDown*pfmetcorr_ex_resolutionDown+pfmetcorr_ey_resolutionDown*pfmetcorr_ey_resolutionDown));
      //if(make_met_plot==true) { cout<<"filling met_resolution_down"<<endl; plotFill("met_resolution_down", mymet_resolutionDown.Pt(), 20, 0, 200,  1.0);}
      return mymet_resolutionDown;
    }
  
}
/* TLorentzVector etau_analyzer::metClusteredUnc(TLorentzVector nominal_met ){ */
/*   TLorentzVector new_metP4; */
/*   TLorentzVector raw_TauP4; */
/*   TLorentzVector uncorrectedMetPlusTau; TLorentzVector met_T1UES; */
/*   TLorentzVector raw_tau = applyTauESCorrections(my_tauP4, TauIndex, 0); */
/*   raw_TauP4.SetPtEtaPhiE(tau_Pt->at(TauIndex),tau_Eta->at(TauIndex), */
/* 			 tau_Phi->at(TauIndex), tau_Energy->at(TauIndex) */
/* 			 );                                                     //// raw tau */
  
/*   if(unc_shift=="up"){ */
/*     met_T1UES.SetPtEtaPhiE(pfMET_T1UESUp ,0, pfMETPhi_T1UESUp, pfMET_T1UESUp); //// raw met */
/*     uncorrectedMetPlusTau = met_T1UES + raw_TauP4; */
/*     TLorentzVector corrected_tau = applyTauESCorrections(raw_TauP4, TauIndex, 0); */
/*     met_T1UES = uncorrectedMetPlusTau - corrected_tau; */
/*     met_T1UES = MetRecoilCorrections(EleIndex, TauIndex, met_T1UES); */
    
/*     double up_percent = abs(met_T1UES.Pt() - nominal_met.Pt())/nominal_met.Pt(); */
/*     if (up_percent>0.2) { up_percent= 0.2; met_T1UES = nominal_met * 1.2;} */
/*     //cout<<"up_percent          "<< up_percent <<endl;    */
/*     new_metP4 = met_T1UES; */
/*   } */
/*   else if(unc_shift=="down") */
/*     { */
/*       met_T1UES.SetPtEtaPhiE(pfMET_T1UESDo ,0, pfMETPhi_T1UESDo, pfMET_T1UESDo); //// raw met */
/*       uncorrectedMetPlusTau = met_T1UES + raw_TauP4; */
/*       TLorentzVector corrected_tau = applyTauESCorrections(raw_TauP4, TauIndex, 0); */
/*       met_T1UES = uncorrectedMetPlusTau - corrected_tau; */
/*       met_T1UES = MetRecoilCorrections(EleIndex, TauIndex, met_T1UES); */

/*       double dn_percent = abs(met_T1UES.Pt() - nominal_met.Pt())/nominal_met.Pt(); */
/*       if (dn_percent>0.2) {dn_percent= 0.2;  met_T1UES = nominal_met * 0.8;} */
/*       new_metP4 = met_T1UES;  */
/*       //cout<<"                                dn_percent "<< dn_percent <<endl; */
/*     } */
/*   else  */
/*     new_metP4 = nominal_met; */
/*   return new_metP4; */
/* } */

TLorentzVector etau_analyzer::metClusteredUnc(TLorentzVector nominal_met ){

  TLorentzVector met_T1UES;

  /// get percent change
  double up_percent = abs(pfMET_T1UESUp - pfMET) / pfMET ;
  double dn_percent = abs(pfMET_T1UESDo - pfMET) / pfMET ;
  double avg_change = (up_percent + dn_percent) / 2;
  if (unc_shift=="up")
    met_T1UES = nominal_met * (1 + avg_change);
  else if (unc_shift=="down")
    met_T1UES = nominal_met * (1 - avg_change);
  else 
    met_T1UES = nominal_met;

  return met_T1UES;
}
