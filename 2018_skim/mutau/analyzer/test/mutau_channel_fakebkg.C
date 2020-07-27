


void mutau_analyzer::Loop_fakebkg(Long64_t maxEvents, int reportEvery, string SampleName, string _isMC_)
{
  
  int nTotal;
  nTotal = 0;
  int report_=0;
  int report_test=0;
  double numberOfEvents = 0;
  int nInspected;
  nInspected = 0;
  double nInspected_genWeighted;  
  nInspected_genWeighted = 0.0; 
  bool debug=false;  
  double netWeight = 1.0;
  double afterSF1=0;
  double afterSF2=0;     
  double afterSF3=0;     
  double afterSF4=0;     

  if (fChain == 0) return;
  //vector<UShort_t> genMatching; genMatching.clear();
  int genMatching=0;
  int thirdLeptonIndex=-1;
  std::vector<int> muonGen;       muonGen.clear();
  std::vector<int> muCand;        muCand.clear();
  std::vector<int> mu2Cand;       mu2Cand.clear();
  std::vector<int> tauCand;       tauCand.clear();
  std::vector<int> jetCand;       jetCand.clear();
  std::vector<int> higgsCand;     higgsCand.clear();
  std::vector<int> reco_mu;       reco_mu.clear(); 
  std::vector<int> reco_mu2;      reco_mu2.clear();
  std::vector<int> reco_tau;      reco_tau.clear(); 
    
  TString sample = TString(SampleName);
  
  int nHiggs = 0;
   int nHToMuTau = 0;
   int found_mt = 0;
   int muCand_1=0; int muCand_2=0;int muCand_3=0;
   int tauCand_1=0; int tauCand_2=0;int tauCand_3=0;

   bool fill_hist = false;
   bool isMC = false;
   if( _isMC_=="MC" ) { isMC=true; fill_hist=true; }
   else if ( _isMC_=="DATA" ) { isMC=false; fill_hist=false; }
  
   Double_t  Pt_Bins[26]={0.0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};
   Double_t  Pt_Bins_highPt[21]={100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};

   //TH1F* h_cutflow=new TH1F("cutflow", "cutflow", 10, 0, 10); h_cutflow->Sumw2();
   TH1F* h_cutflow_n=new TH1F("cutflow_n", "cutflow_n", 8, 0, 8);h_cutflow_n->Sumw2();
   TH1F* h_cutflow_n_fr=new TH1F("cutflow_n_fr", "cutflow_n_fr", 8, 0, 8);h_cutflow_n_fr->Sumw2();
   TH1F* h_cutflow_n_dyll=new TH1F("cutflow_n_dyll", "cutflow_n_dyll", 8, 0, 8);h_cutflow_n_dyll->Sumw2();
   TH1F* h_cutflow_n_dyll_fr=new TH1F("cutflow_n_dyll_fr", "cutflow_n_dyll_fr", 8, 0, 8);h_cutflow_n_dyll_fr->Sumw2();
   //TH1F* h_cutflow_Htt=new TH1F("cutflow_Htt", "cutflow_Htt", 11, 0, 11); h_cutflow_Htt->Sumw2();

   TLorentzVector myMomTau, myTauh,  myNeu, myHiggs; 
   bool found_Wjet_sample=false;
   bool found_DYjet_sample=false;
   int PID =0;
   if ( sample.Contains("WJetsToLNu") ||
	sample.Contains("W1JetsToLNu") ||
	sample.Contains("W2JetsToLNu") ||
	sample.Contains("W3JetsToLNu") ||
	sample.Contains("W4JetsToLNu") 	) {
     found_Wjet_sample=true;
     PID = 24;
     cout<<"****************** wjet sample found"<<endl;
   }
   if ( sample.Contains("DYJetsToLL") ||
	sample.Contains("DY1JetsToLL") ||
	sample.Contains("DY2JetsToLL") ||
	sample.Contains("DY3JetsToLL") ||
	sample.Contains("DY4JetsToLL")  ) {
     found_DYjet_sample=true;
     PID = 23;
     cout<<"****************** dyjet sample found"<<endl;
   }
   if(debug)cout<<" setting up kFactor files ..."<<endl;
   TH1F *NLO_QCD_EWK,*NLO_EWK,*NLO_QCD,*NNLO_QCD;
   TFile* f_nnlo_qcd = TFile::Open("sf_files/RootFiles/theory/lindert_qcd_nnlo_sf.root");
   TFile* f_nlo_qcd  = TFile::Open("sf_files/RootFiles/theory/2017_gen_v_pt_qcd_sf.root");
   TFile* f_qcd_ewk;
   if ( found_Wjet_sample ) {
     f_qcd_ewk = TFile::Open("sf_files/RootFiles/theory/merged_kfactors_wjets.root");
     NLO_QCD_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_qcd_ewk");
     NLO_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_ewk");
     NLO_QCD = (TH1F*)f_nlo_qcd->Get("wjet_dress_monojet");
     NNLO_QCD = (TH1F*)f_nnlo_qcd->Get("evj");
     
   } else if ( found_DYjet_sample ) {
     f_qcd_ewk = TFile::Open("sf_files/RootFiles/theory/merged_kfactors_zjets.root");
     NLO_QCD_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_qcd_ewk");
     NLO_EWK = (TH1F*)f_qcd_ewk->Get("kfactor_monojet_ewk");
     f_nlo_qcd = TFile::Open("sf_files/RootFiles/theory/kfac_dy_filter.root");
     NLO_QCD = (TH1F*)f_nlo_qcd->Get("kfac_dy_filter");
     NNLO_QCD = (TH1F*)f_nnlo_qcd->Get("eej");
     
   }
   
   // if(debug)cout<<" setting up other files ..."<<endl;
   RoccoR  rc("sf_files/roCorr_Run2_v3/RoccoR2017.txt"); 
   TLorentzVector muonP4;
   TLorentzVector tauP4;
   
   Long64_t nentries = fChain->GetEntries();
   if ( isMC==true ) std::cout<<".... MC file ..... "<<std::endl;
   else  std::cout<<".... DATA file ..... "<<std::endl;

   std::cout<<"Coming in: "<<std::endl;
   std::cout<<"nentries:"<<nentries<<std::endl;
   //Look at up to maxEvents events, or all if maxEvents == -1.
   Long64_t nentriesToCheck = nentries;
   if (maxEvents != -1LL && nentries > maxEvents)
     nentriesToCheck = maxEvents;
   nTotal = nentriesToCheck;
   Long64_t nbytes = 0, nb = 0;
   
   std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
   //TStopwatch sw;
   //sw.Start();
   
   for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
     {
       if(debug) cout<<"event "<<jentry<<endl;
       muCand.clear(); mu2Cand.clear();
       tauCand.clear();
       jetCand.clear();
       reco_mu.clear(); 
       reco_mu2.clear();
       reco_tau.clear();  
       muonGen.clear();
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       double inspected_event_weight = 1.0; 
       if(isMC)	 fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight/fabs(genWeight) : inspected_event_weight = 0.0;
       nInspected_genWeighted += inspected_event_weight;  
       nInspected += 1; 
       //h_insEvents->SetBinContent(1, nInspected_genWeighted);
       //=1.0 for real data
       double event_weight=1.0;
       double weight=1.0;
       double muRC_sf = 1.0; double randomN = gRandom->Rndm();
       double sf_tauID = 1.0; 

       double pileup_sf = 1.0;
       double kfactor = 1.0;
       double nlo_ewk = 1.0;
       double nlo_qcd_binned = 1.0;
       double nlo_qcd = 1.0;
       double nnlo_qcd =1.0;
       int bosonPID;
       double bosonPt=0.0;
       bool Wfound = false;
       bool passSingleTriggerPaths=false;
       bool passCrossTrigger=false;
       int report_i=0;
       bool muTriggerFilterMatch=false;
       bool tauTriggerFilterMatch=false;
       bool muTau_selector=false;
       bool muMu_selector=false;
       numberOfEvents+=weight;
       weight=inspected_event_weight;
       
       if(debug)cout<<"this worked Line 314"<<endl;
       
       if(isMC) weight=inspected_event_weight;
       else weight=1.0;
       int leading_muon = -1; float leading_mPt=0;
       int leading_tau = -1;  float leading_tPt=0;
       
       if(isMC)
	 pileup_sf = h_pileup->GetBinContent(h_pileup->GetXaxis()->FindBin(puTrue->at(0)));
       weight = weight*pileup_sf;
       if( isGoodVtx==false ) continue;
       //if( noisyJet2017()==true ) continue;
       
       if( found_DYjet_sample && hasGenTau())
	 muTau_selector=true;
       else if( found_DYjet_sample && !hasGenTau() )
	 muTau_selector=false;
       else if ( !found_DYjet_sample )
	 muTau_selector=true;

       if( found_DYjet_sample && 
	   (!found_GenMatch(5) && !found_GenMatch(6))
	   )
	 muMu_selector=true;
       else
	 muMu_selector=false;
       
       /////Trigger bit selection
       if(HLTEleMuX>>21&1 == 1 || HLTEleMuX>>60&1 == 1 )
	 passSingleTriggerPaths=true;
       if(  (HLTTau>>13&1 == 1 && !isMC)
	    || HLTTau>>14&1 == 1
	    )
	 passCrossTrigger=true;
       ////
       if(isMC && (found_Wjet_sample || found_DYjet_sample)){
	 if(debug)cout<<"check which mc particle is W boson"<<endl;
	 for(int i=0; i<nMC;i++){
	   if(abs((*mcPID)[i]) == PID){
	     if(!( mcStatus->at(i) == 62)) continue;
	     Wfound=true;
	     bosonPID = (*mcPID)[i];
	     bosonPt = (*mcPt)[i];
	   }
	 }
	 if ( bosonPt > 0 ){
	   if(debug)cout<<"Accessing nlo ewk, qcd"<<endl;
	   nlo_ewk = NLO_EWK->GetBinContent(NLO_EWK->GetXaxis()->FindBin(bosonPt));
	   nlo_qcd_binned=NLO_QCD->GetBinContent(NLO_QCD->GetXaxis()->FindBin(bosonPt));
	   if (found_Wjet_sample) {
	     nlo_qcd = exponential(bosonPt,1.053, 3.163e-3, 0.746);
	   } else if (found_DYjet_sample) {
	     //nlo_qcd = exponential(bosonPt,1.434, 2.210e-3, 0.443);
	     nlo_qcd = NLO_QCD->GetBinContent(NLO_QCD->GetXaxis()->FindBin(bosonPt));
	   } //else if (type == GJets) {
	   //nlo_qcd = exponential(bosonPt,1.159, 1.944e-3, 1.0);
	   // }
	   
	   if(debug)cout<<"Accessing nnlo qcd"<<endl;
	   nnlo_qcd = NNLO_QCD->GetBinContent(NNLO_QCD->GetXaxis()->FindBin(bosonPt));
	 }
	 // if (isNLO) kfactor = nlo_ewk * nnlo_qcd;
	 // else kfactor = nlo_ewk * nlo_qcd * nnlo_qcd;
      	 if(debug) cout<<"apply kfactor"<<endl;
	 if(nlo_ewk*nlo_qcd !=0 ) kfactor = nlo_ewk * nlo_qcd;
	 
       }
       weight=weight*kfactor;
       //cout<<" kfactor = "<<kfactor<<endl;
       event_weight=weight;
       if(debug)cout<<"reco selections begin"<<endl;
       
	 
       ////// fake background region - antiisolated begin
       if(debug)cout<<"moving to fake bkg"<<endl;
       event_weight=weight;
       muCand.clear(); tauCand.clear();
       if(metFilters==0)
	 {
	   
	   if(isMC) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
	   nMETFiltersPassed_fr+=event_weight;
	   if(  passSingleTriggerPaths || passCrossTrigger  )
	     {
	       nSingleTrgPassed_fr+=event_weight;
	       if(debug)cout<<"trigger selected"<<endl;
	       muCand = getMuCand(20,2.1);  ///// muons selected 
	       if( muCand.size() >0 ) 
		 { 
		   nGoodMuonPassed_fr+=event_weight;
		   if(debug)cout<<"this worked Line 373"<<endl;
		   tauCand = getAISRTauCand(30,2.3);
		   if( tauCand.size()>0 ) 
		     {
		       nGoodTauPassed_fr+=event_weight;
		       if(debug)cout<<"fr tau selection passed"<<endl;
		       
		       reco_mu.clear();reco_tau.clear();
		       reco_mu=muCand; reco_tau=tauCand;
		       if ( MatchTriggerFilter(reco_mu[0], reco_tau[0]) )
			 {
			   
			   if ( muCharge->at(reco_mu[0]) * tau_Charge->at(reco_tau[0]) < 0  ) 
			     {
			       nGoodMuTauPassed_fr+=event_weight;
			       muonP4.SetPtEtaPhiE(muPt->at(reco_mu[0]), muEta->at(reco_mu[0]), muPhi->at(reco_mu[0]), muE->at(reco_mu[0]));
			       tauP4.SetPtEtaPhiE(tau_Pt->at(reco_tau[0]), tau_Eta->at(reco_tau[0]), tau_Phi->at(reco_tau[0]), tau_Energy->at(reco_tau[0]));
			       
			       
			       muonP4 = muonP4*muRC_sf;
			       
			       event_weight = event_weight* getFR(reco_tau[0]);
			       
			       /////
			       if(debug)cout<<" fake bkg sf : "<<getScaleFactors( reco_mu[0] , reco_tau[0] , true , isMC , debug ) <<endl;
			       if(isMC) event_weight = event_weight * getScaleFactors( reco_mu[0] , reco_tau[0] , true , isMC , debug );
			       /////
			       if( thirdLeptonVeto() < 0 )
				 {
				   nPassedThirdLepVeto_fr+=event_weight;
				   
				   if( passBjetVeto() == true)
				     {
				       nPassedBjetVeto_fr+=event_weight;
				       
				       double deltaR = delta_R(muonP4.Phi(), muonP4.Eta(), tauP4.Phi(), tauP4.Eta());
				       if(deltaR > 0.5 )
					 {
					   nDeltaRPassed_fr+=event_weight;
					   if(debug)cout<<"this worked Line 425"<<endl;
					   if(debug)cout<<"evnt_weight = "<<event_weight<<endl;
					   fillHist("5_fr", muonP4, tauP4, reco_mu[0], reco_tau[0], event_weight);
					   double mT_muMet = TMass_F((muPt->at(reco_mu[0])),(muPhi->at(reco_mu[0])),pfMET,pfMETPhi  );
					   if( mT_muMet < 50 ) 
					     {
					       fillHist("6_fr", muonP4, tauP4, reco_mu[0], reco_tau[0], event_weight);
					     }
					 }
				     }
				 }
			     }
			 }
		     }
		 }
	     }
	 }
     
       ////// fake rate anti isolated region end

       
       report_test = nentriesToCheck/20;
       while (report_test>10)
	 {
	   report_test=report_test/10;
	   report_i++;
	 }
       if(nentriesToCheck>20)
	 reportEvery = report_test*pow(10,report_i);
       else 
	 reportEvery = 1;
       if (jentry%reportEvery == 0)
	 {
	   std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
	 }
	   
     }
   

   
   h_cutflow_n_fr->SetBinContent(1,nInspected_genWeighted );
   h_cutflow_n_fr->SetBinContent(2, nSingleTrgPassed_fr);
   h_cutflow_n_fr->SetBinContent(3, nGoodMuonPassed_fr);
   h_cutflow_n_fr->SetBinContent(4, nGoodTauPassed_fr);
   h_cutflow_n_fr->SetBinContent(5, nGoodMuTauPassed_fr);
   h_cutflow_n_fr->SetBinContent(6, nPassedThirdLepVeto_fr);
   h_cutflow_n_fr->SetBinContent(7, nPassedBjetVeto_fr);
   h_cutflow_n_fr->SetBinContent(8, nDeltaRPassed_fr);

   fileName->cd();
   map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
   map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
   for (; iMap1 != jMap1; ++iMap1)
     nplot1(iMap1->first)->Write();
}
