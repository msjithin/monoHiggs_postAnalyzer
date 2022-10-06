#ifndef HTT_MEtSys_h
#define HTT_MEtSys_h

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TString.h>
#include <TRandom.h>
#include <TMath.h>
#include <assert.h>

class MEtSys {
  
 public:
  MEtSys(){
    /* TString cmsswBase = TString( getenv ("CMSSW_BASE") ); */
    /* TString baseDir = cmsswBase + "/src"; */
    /* TString _fileName = baseDir+"/"+fileName; */
    TString fileName = "sf_files/HTT-utilities/RecoilCorrections/data/PFMEtSys_2017.root";
    TString _fileName = "sf_files/HTT-utilities/RecoilCorrections/data/PFMEtSys_2017.root";
    TFile * file = new TFile(_fileName);
    if (file->IsZombie()) {
      std::cout << "file " << _fileName << " is not found...   quitting " << std::endl;
      exit(-1);
    }
    TH1D * typeBkgdH = (TH1D*)file->Get("typeBkgdH");
    if (typeBkgdH==NULL) {
      std::cout << "Histogram typeBkgdH should be contained in file " << fileName << std::endl;
      std::cout << "Check content of the file " << _fileName << std::endl;
      exit(-1);
    }
    TH1D * jetBinsH = (TH1D*)file->Get("jetBinsH");
    if (jetBinsH==NULL) {
      std::cout << "Histogram jetBinsH should be contained in file " << fileName << std::endl;
      std::cout << "Check content of the file " << _fileName << std::endl;
      exit(-1);
    }

    nBkgdTypes = typeBkgdH->GetNbinsX();
    std::vector<TString> Bkgd; Bkgd.clear();
    nJetBins = jetBinsH->GetNbinsX();
    std::vector<TString> JetBins; JetBins.clear();
    for (int i=0; i<nBkgdTypes; ++i) {
      Bkgd.push_back(typeBkgdH->GetXaxis()->GetBinLabel(i+1));
    }
    for (int i=0; i<nJetBins; ++i) {
      JetBins.push_back(jetBinsH->GetXaxis()->GetBinLabel(i+1));
    }

    TString uncType[2] = {"Response",
			  "Resolution"};

    for (int i=0; i<nBkgdTypes; ++i) {
      TString histName = Bkgd[i]+"_syst";
      TH2D * hist = (TH2D*)file->Get(histName);
      if (hist==NULL) {
	std::cout << "Histogram " << histName << " should be contained in file " << fileName << std::endl;
	std::cout << "Check content of the file " << fileName << std::endl;
	exit(-1);
      }
      for (int xBin=0; xBin<2; ++xBin) {
	for (int yBin=0; yBin<3; ++yBin) {
	  sysUnc[i][xBin][yBin] = hist->GetBinContent(xBin+1,yBin+1);
	  std::cout << "Systematics : " << Bkgd[i] << "  " << uncType[xBin] << " " << JetBins[yBin] << " = " << sysUnc[i][xBin][yBin] << std::endl;
	}
      }

    }
    for (int i=0; i<nBkgdTypes; ++i) {
      for (int j=0; j<nJetBins; ++j) {
	TString histName = Bkgd[i]+"_"+JetBins[j];
	std::cout << histName << std::endl;
	responseHist[i][j] = (TH1D*)file->Get(histName);
	if (responseHist[i][j]==NULL) {
	  std::cout << "Histogram " << histName << " should be contained in file " << fileName << std::endl;
	  std::cout << "Check content of the file " << fileName << std::endl;
	  exit(-1);
	}

      }
    }
  }

  ~MEtSys(){};

  void ApplyMEtSys(float metPx,
		   float metPy,
		   float genVPx,
		   float genVPy,
		   float visVPx,
		   float visVPy,
		   int njets,
		   int bkgdType,
		   int sysType,
		   int sysShift,
		   float & metShiftPx,
		   float & metShiftPy) {

    if (bkgdType<0||bkgdType>=nBkgdTypes) {
      std::cout << "MEtSys::ShiftResponseMet() : Background type " << bkgdType << " does not correspond to any of allowed options : " << std::endl;
    std::cout << "0 : Z(W)+Jets" << std::endl;
    std::cout << "1 : EWK+single-top" << std::endl;
    std::cout << "2 : top pair" << std::endl;
    exit(-1);
  }

  int jets = njets;
  if (jets>2) jets = 2;
  if (jets<0) {
    std::cout << "MEtSys::ApplyMEtSys() : Number of jets is negative !" << std::endl;
    exit(-1);
  }

  int type = 0; if (sysType!=0) type = 1;

  float scale = 1 + sysUnc[bkgdType][type][jets];
  if (sysShift!=0) scale = 1 - sysUnc[bkgdType][type][jets];

  ShiftMEt(metPx,
           metPy,
           genVPx,
           genVPy,
           visVPx,
           visVPy,
           njets,
           bkgdType,
           type,
           scale,
           metShiftPx,
           metShiftPy);
}

  void ShiftMEt(float metPx,
float metPy,
float genVPx, 
float genVPy,
float visVPx,
float visVPy,
int njets,
int bkgdType,
int sysType,
float sysShift,
float & metShiftPx,
float & metShiftPy){
  metShiftPx=metPx;
  metShiftPy=metPy;

  if (sysType==0)
    ShiftResponseMet(metPx,
                     metPy,
                     genVPx,
                     genVPy,
                     visVPx,
                     visVPy,
                     njets,
                     bkgdType,
                     sysShift,
                     metShiftPx,
                     metShiftPy);
  else
    ShiftResolutionMet(metPx,
                       metPy,
                       genVPx,
                       genVPy,
                       visVPx,
                       visVPy,
                       njets,
                       bkgdType,
                       sysShift,
                       metShiftPx,
                       metShiftPy);
}

  void ShiftResponseMet(float metPx,
float metPy,
float genVPx, 
float genVPy,
float visVPx,
float visVPy,
int njets,
int bkgdType,
float sysShift,
float & metShiftPx,
float & metShiftPy){
  float Hparal = 0;
  float Hperp = 0;
  float genVPt = TMath::Sqrt(genVPx*genVPx+genVPy*genVPy);

  // protection against null
  if (genVPt<1.0) {
    metShiftPx = metPx;
    metShiftPy = metPy;
    return;
  }

  ComputeHadRecoilFromMet(metPx,metPy,genVPx,genVPy,visVPx,visVPy,Hparal,Hperp);

  int jets = njets;
  if (jets>2) jets = 2;
  if (jets<0) {
    std::cout << "MEtSys::ShiftResponseMet() : Number of jets is negative !" << std::endl;
    exit(-1);
  }

  if (bkgdType<0||bkgdType>=nBkgdTypes) {
    std::cout << "MEtSys::ShiftResponseMet() : Background type " << bkgdType << " does not correspond to any of allowed options : " << std::endl;
    std::cout << "0 : Z(W)+Jets" << std::endl;
    std::cout << "1 : EWK+single-top" << std::endl;
    std::cout << "2 : top pair" << std::endl;
    exit(-1);
  }
  float mean = -responseHist[bkgdType][jets]->Interpolate(genVPt)*genVPt;
  float shift = sysShift*mean;
  Hparal = Hparal + (shift-mean);

  ComputeMetFromHadRecoil(Hparal,Hperp,genVPx,genVPy,visVPx,visVPy,metShiftPx,metShiftPy);

  }

  
  void ShiftResolutionMet(float metPx,
			  float metPy,
			  float genVPx, 
			  float genVPy,
			  float visVPx,
			  float visVPy,
			  int njets,
			  int bkgdType,
			  float sysShift,
			  float & metShiftPx,
			  float & metShiftPy){
    float Hparal = 0;
    float Hperp = 0;
    float genVPt = TMath::Sqrt(genVPx*genVPx+genVPy*genVPy);

    // protection against null
    if (genVPt<1.0) {
      metShiftPx = metPx;
      metShiftPy = metPy;
      return;
    }

    ComputeHadRecoilFromMet(metPx,metPy,genVPx,genVPy,visVPx,visVPy,Hparal,Hperp);

    int jets = njets;
    if (jets>2) jets = 2;
    if (jets<0) {
      std::cout << "MEtSys::ShiftResponseMet() : Number of jets is negative !" << std::endl;
      exit(-1);
    }

    if (bkgdType<0||bkgdType>=nBkgdTypes) {
      std::cout << "MEtSys::ShiftResponseMet() : Background type " << bkgdType << " does not correspond to any of allowed options : " << std::endl;
      std::cout << "0 : Z(W)+Jets" << std::endl;
      std::cout << "1 : EWK+single-top" << std::endl;
      std::cout << "2 : top pair" << std::endl;
      exit(-1);
    }
    float mean = -responseHist[bkgdType][jets]->Interpolate(genVPt)*genVPt;
    Hperp = sysShift*Hperp;
    Hparal = mean + (Hparal-mean)*sysShift;

    ComputeMetFromHadRecoil(Hparal,Hperp,genVPx,genVPy,visVPx,visVPy,metShiftPx,metShiftPy);

  }

  enum ProcessType{BOSON=0, EWK=1, TOP=2};
  enum SysType{Response=0, Resolution=1};
  enum SysShift{Up=0, Down=1};

 private:

  void ComputeHadRecoilFromMet(float metX,
			       float metY,
			       float genVPx, 
			       float genVPy,
			       float visVPx,
			       float visVPy,
			       float & Hparal,
			       float & Hperp){
    float genVPt = TMath::Sqrt(genVPx*genVPx+genVPy*genVPy);
    float unitX = genVPx/genVPt;
    float unitY = genVPy/genVPt;

    float unitPhi = TMath::ATan2(unitY,unitX);
    float unitPerpX = TMath::Cos(unitPhi+0.5*TMath::Pi());
    float unitPerpY = TMath::Sin(unitPhi+0.5*TMath::Pi());

    float Hx = -metX - visVPx;
    float Hy = -metY - visVPy;

    Hparal = Hx*unitX + Hy*unitY;
    Hperp = Hx*unitPerpX + Hy*unitPerpY;
  }


  void ComputeMetFromHadRecoil(float Hparal,
			       float Hperp,
			       float genVPx, 
			       float genVPy,
			       float visVPx,
			       float visVPy,
			       float & metX,
			       float & metY){
    float genVPt = TMath::Sqrt(genVPx*genVPx+genVPy*genVPy);
    float unitX = genVPx/genVPt;
    float unitY = genVPy/genVPt;

    float unitPhi = TMath::ATan2(unitY,unitX);
    float unitPerpX = TMath::Cos(unitPhi+0.5*TMath::Pi());
    float unitPerpY = TMath::Sin(unitPhi+0.5*TMath::Pi());

    float det = unitX*unitPerpY - unitY*unitPerpX;
    float Hx = (Hparal*unitPerpY - Hperp*unitY)/det;
    float Hy = (Hperp*unitX - Hparal*unitPerpX)/det;

    metX = -Hx - visVPx;
    metY = -Hy - visVPy;
  }
  
  
  int nBkgdTypes;
  int nJetBins;
  TH1D * responseHist[3][5];
  float sysUnc[3][2][3]; // first  index : bkgd type 
  // second index : type of uncertainty 0=response, 1=resolution
  // third index  : jet multiplicity bin (0,1,2);

};


#endif
