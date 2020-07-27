void etau_analyzer::BookHistos(const char* file1, const char* file2)
{
  TFile* file_in= new TFile(file1, "READ");
  fileName = new TFile(file2, "RECREATE");
  
  //makeOutputTree(tree);
  fileName->cd();
  h_nEvents = (TH1F*)((TH1F*)file_in->Get("nEvents"))->Clone(TString("nEvents"));
  file_in->Close();
  

}

//Fill the sequential histos at a particular spot in the sequence


void etau_analyzer::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{
  
  //  h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);

}

double etau_analyzer::dR(int ele_index, int tau_index)
{
  double deltaeta = abs(eleEta->at(ele_index) - tau_Eta->at(tau_index));
  double electronPhi = elePhi->at(ele_index);
  double tauPhi = tau_Phi->at(tau_index);

  double deltaphi = DeltaPhi(electronPhi, tauPhi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}

double etau_analyzer::delta_R(float phi1, float eta1, float phi2, float eta2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar   = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
  
}



double etau_analyzer::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
}

float etau_analyzer::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
  return  sqrt(2.0*LepPt*met*(1.0-cos(DeltaPhi(LepPhi, metPhi))));
  //return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

float etau_analyzer::TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float totalTMass = (a + b+ met).M();
  return totalTMass;
}


float etau_analyzer::VisMass_F(TLorentzVector a, TLorentzVector b){
  float visibleMass = (a + b).M();
  return visibleMass;
}

float etau_analyzer::pTvecsum_F(float pt1, float pt2, float phi1, float phi2) {
  float pt_vecSum = sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2), 2) + pow(pt1*sin(phi1) + pt2*sin(phi2), 2));
  return pt_vecSum;
}
float etau_analyzer::pTvecsum_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float pt_vecSum = (a + b+ met).Pt();
  return pt_vecSum;
}

bool etau_analyzer::passBjetVeto()
{
  std::vector<int> tmpJetCand;
  tmpJetCand.clear();
  bool veto = true;
  bool foundBjet = false;
  // for(int iJets=0; iJets<nJet ; iJets++){
  //   if(jetCSV2BJetTags->at(iJets) > 0.8838) tmpCand++;
  //   // CSV B jet tag for selecting bJets is medium WP (jetCSV2BJetTags > 0.8838.)
  // }
  for(int iJets=0; iJets<nJet ; iJets++){
    if( jetPt->at(iJets) > 25  && abs(jetEta->at(iJets)) < 2.4 && jetID->at(iJets)>>1&1==1 ){
      tmpJetCand.push_back(iJets);
    }
  }
  if(tmpJetCand.size() >=1 ){
    // atleast one jet ==> events pass medium 
    if( (jetDeepCSVTags_b->at(tmpJetCand[0]) + jetDeepCSVTags_bb->at(tmpJetCand[0])) > 0.4941  )
      return veto = false;
  }
  else if(tmpJetCand.size() >= 2){
    // atleast 2 jets ==> events pass loose
    if( (jetDeepCSVTags_b->at(tmpJetCand[0]) + jetDeepCSVTags_bb->at(tmpJetCand[0])) > 0.1522   )
      return veto = false;
  }
  return veto;
}
bool etau_analyzer::passDiElectronVeto(int eleIndex)
{
  std::vector<int> tmpCand; int tmpEleIndex1=-1; int tmpEleIndex2=-1;
  tmpCand.clear();
  bool veto = true;
  bool awayFromEverything=true;
  for(int iEle=0;iEle<nEle;iEle++)
    {
      bool kinematic = false;
      if( elePt->at(iEle) > 15
	  && fabs(eleEta->at(iEle))< 2.5
	  && fabs(eleD0->at(iEle)) < 0.045
	  && fabs(eleDz->at(iEle)) < 0.2
       	  ) kinematic = true;
      bool electronId =false;
      if( eleIDbit->at(iEle)>>3&1==1) electronId =true; // cut based electron id veto
      bool relative_iso = false;    
      float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
      if( relEleIso < 0.3 ) relative_iso = true;
      if( kinematic && electronId && relative_iso ){
	tmpCand.push_back(iEle);
      }	
    }
  std::vector<int> iElePlus; iElePlus.clear(); 
  std::vector<int> iEleMinus; iEleMinus.clear();
  for(int i=0; i<tmpCand.size(); i++){
    if(eleCharge->at(tmpCand[i]) < 0) iEleMinus.push_back(tmpCand[i]);
    if(eleCharge->at(tmpCand[i]) > 0) iElePlus.push_back(tmpCand[i]);
  }
  
  if(iElePlus.size()>0 && iEleMinus.size()>0){
    for(int i=0;i<iElePlus.size();i++)
      {
	for(int j=0;j<iEleMinus.size();j++)
	  {
	    double deltaR= delta_R(elePhi->at(iEleMinus[j]), eleEta->at(iEleMinus[j]), elePhi->at(iElePlus[i]), eleEta->at(iElePlus[i]));
	    if (deltaR < 0.15) {
	      awayFromEverything = false; tmpEleIndex1=i; tmpEleIndex2=j;
	      break;
	    }
	  }
      }
    if (awayFromEverything  && eleCharge->at(iElePlus[0])*eleCharge->at(iEleMinus[0])<0) {
      return false;
    }
  }
  
  return true;
  
}
// bool etau_analyzer::passDiElectronVeto(int eleIndex)
// {
//   std::vector<int> tmpCand; int tmpEleIndex=-1;
//   tmpCand.clear();
//   bool veto = true;
//   bool awayFromEverything=true;
//   for(int iEle=0;iEle<nEle;iEle++)
//     {
//       bool kinematic = false;
//       if( elePt->at(iEle) > 15
// 	  && fabs(eleEta->at(iEle))< 2.5
// 	  && fabs(eleD0->at(iEle)) < 0.045
// 	  && fabs(eleDz->at(iEle)) < 0.2
//        	  ) kinematic = true;
//       bool electronId =false;
//       if( eleIDbit->at(iEle)>>3&1==1) electronId =true; // cut based electron id veto
//       bool relative_iso = false;    
//       float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
//       if( relEleIso < 0.3 ) relative_iso = true;
//       if( kinematic && electronId && relative_iso ){
// 	tmpCand.push_back(iEle);
//       }	
//     }
//   if(tmpCand.size()>0){
//     for(int iEle=0;iEle<tmpCand.size();iEle++)
//       {
// 	awayFromEverything=true;
// 	double deltaR= delta_R(elePhi->at(eleIndex), eleEta->at(eleIndex), elePhi->at(tmpCand[iEle]), eleEta->at(tmpCand[iEle]));
// 	if (deltaR < 0.15) {
// 	  awayFromEverything = false; tmpEleIndex=iEle;
// 	  break;
// 	}
//       }
//     if (!awayFromEverything && tmpEleIndex>-1 && eleCharge->at(eleIndex)*eleCharge->at(tmpEleIndex)<0) {
//       return false;
//     }
//   }
  
//   return true;
  
// }
int etau_analyzer::eVetoZTTp001dxyz(double minDeltaR){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool awayFromEverything = true;   int tmpEleIndex=-1;
  //Loop over electrons      
  for(int iEle=0;iEle<nEle;iEle++)
    {
      bool kinematic = false;
      if( elePt->at(iEle) > 10
	  && fabs(eleEta->at(iEle))< 2.5
	  && fabs(eleD0->at(iEle)) < 0.045
	  && fabs(eleDz->at(iEle)) < 0.2
	  && eleConvVeto->at(iEle)==1
       	  ) kinematic = true;
      bool electronId =false;
      if( eleIDbit->at(iEle)>>8&1==1) electronId =true; // cut based electron id veto
      bool relative_iso = false;    
      float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
      if( relEleIso < 0.3 ) relative_iso = true;
      if( kinematic && electronId && relative_iso ){
	tmpCand.push_back(iEle);
      }	
    }
  std::vector<int> output;         output.clear();
  std::vector<int> tauCand;        tauCand.clear();
  tauCand = getTauCand(30,2.3);
  for (int i = 0; i < tmpCand.size(); ++i) {
    bool awayFromEverything = true;
    for (int j = 0; j < tauCand.size(); ++j) {
      double deltaR = delta_R(tau_Phi->at(tauCand[j]), tau_Eta->at(tauCand[j]), elePhi->at(tmpCand[i]), eleEta->at(tmpCand[i]));
      if (deltaR < minDeltaR) {
        awayFromEverything = false; tmpEleIndex=tmpCand[i];
        break;
      }
    }
    if (awayFromEverything && tmpCand.size()>0) {
      output.push_back(i);   
    }
  }
  return output.size();
   
}
int etau_analyzer::mVetoZTTp001dxyz(double minDeltaR){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool awayFromEverything = true;   int tmpMuIndex=-1;
  //Loop over muons
  for(int iMu=0; iMu < nMu;iMu++)
    {
      bool kinematic = false;
      if( (*muPt)[iMu] > 10.0  && fabs((*muEta)[iMu])< 2.4 && (*muD0)[iMu] < 0.045 && (*muDz)[iMu] < 0.2 ) kinematic = true;
      bool muonId =false;
      if( muIDbit->at(iMu)>>1&1==1) muonId =true;
      bool relative_iso = false;
      float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
      if( relMuIso < 0.3 ) relative_iso = true;
      if(muonId==true && kinematic==true && relative_iso==true){
	tmpCand.push_back(iMu);
      }                   
    }          

  std::vector<int> output;         output.clear();
  std::vector<int> tauCand;        tauCand.clear();
  tauCand = getTauCand(30,2.3);
  for (int i = 0; i < tmpCand.size(); ++i) {
    bool awayFromEverything = true;
    for (int j = 0; j < tauCand.size(); ++j) {
      double deltaR = delta_R(tau_Phi->at(tauCand[j]), tau_Eta->at(tauCand[j]), muPhi->at(tmpCand[i]), muEta->at(tmpCand[i]));
      if (deltaR < minDeltaR) {
        awayFromEverything = false; tmpMuIndex=tmpCand[i];
        break;
      }
    }
    if (awayFromEverything && tmpCand.size()>0) {
      output.push_back(i);   
    }
  }
  return output.size();
   
}

std::vector<int> etau_analyzer::gen_matching(){
  int tmpCand=-1;
  std::vector<int> tmpGenMatch;
  tmpGenMatch.clear();
  
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>1&1==1 ) tmpGenMatch.push_back(1);
    if( genMatch2->at(imc)>>2&1==1 ) tmpGenMatch.push_back(2);
    if( genMatch2->at(imc)>>3&1==1 ) tmpGenMatch.push_back(3);
    if( genMatch2->at(imc)>>4&1==1 ) tmpGenMatch.push_back(4);
    if( genMatch2->at(imc)>>5&1==1 ) tmpGenMatch.push_back(5);
    if( genMatch2->at(imc)>>6&1==1 ) tmpGenMatch.push_back(6);
  }
  
  // if(tmpGenMatch.size() >0 )
  //   tmpCand=tmpGenMatch[0];
  // return tmpCand; 
  return tmpGenMatch;
}
bool  etau_analyzer::found_GenMatch(int genTau)
{
  std::vector<int> v = gen_matching();
  if (std::find(v.begin(), v.end(), genTau) != v.end())
    return true;
  
  return false;
}
std::vector<int> etau_analyzer::getGenMu(){
  std::vector<int> tmpCand;
  tmpCand.clear();
  int count1=0; int count2=0;
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>2&1==1 || genMatch2->at(imc)>>4&1==1 ) {tmpCand.push_back(imc); count1++;}
    if( (genMatch2->at(imc)>>2&1==1 || genMatch2->at(imc)>>4&1==1) && fabs(mcPID->at(imc))==13 ){count2++;}
  }
  //cout<<"count1:"<<count1<<"  count2:"<<count2<<endl;
  return tmpCand; 
}
bool etau_analyzer::hasGenTau(){
  bool found_genTau=false;
  for(int imc=0; imc<nMC; imc++){
    if( genMatch2->at(imc)>>5&1==1) {  found_genTau=true;}
  }
  return found_genTau;
}
float etau_analyzer::exponential(float x,float a,float b,float c) {
  return a * TMath::Exp(-b * x) + c;
}

double  etau_analyzer::getFR(int tauIndex){
    double frWeight=1.0;
  double tau_FR = 1.0;
  double tauPt=0.0;
  if( tau_Pt->at(tauIndex) < 120 )
    tauPt=tau_Pt->at(tauIndex);
  else
    tauPt=119.0;
  if ( tau_DecayMode->at(tauIndex)==0 )
    {
      tau_FR = h_tauFR_0->Eval(tauPt);
      frWeight = tau_FR/(1-tau_FR);
    }
  
  if ( tau_DecayMode->at(tauIndex)==1 )
    {
      tau_FR = h_tauFR_1->Eval(tauPt);
      frWeight = tau_FR/(1-tau_FR);
    }
  
  if ( tau_DecayMode->at(tauIndex)==10 )
    {
      tau_FR = h_tauFR_10->Eval(tauPt);
      frWeight = tau_FR/(1-tau_FR);
    }
  if ( tau_DecayMode->at(tauIndex)==11 )
    {
      tau_FR = h_tauFR_11->Eval(tauPt);
      frWeight = tau_FR/(1-tau_FR);
    }
  return frWeight;
}
