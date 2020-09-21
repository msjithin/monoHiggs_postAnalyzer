Set up
------

```
mkdir ComputeFF
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git clone https://github.com/cecilecaillol/ComputeFF2018.git
git clone https://github.com/CMS-HTT/LeptonEfficiencies.git
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs
mkdir TauAnalysisTools
git clone -b run2_SFs https://github.com/cms-tau-pog/TauTriggerSFs $CMSSW_BASE/src/TauAnalysisTools/TauTriggerSFs
git clone https://github.com/cms-tau-pog/TauTriggerSFs.git 
scram b -j 8
cd ComputeFF2018
```

To run everything
-----------------

```
cd ComputeFF2018/
MakeFFs.py --year <year, 2016,2017,2018> --channel <mt,et,tt>
```

To change the ID WP:
-------------------

Edit FFcode/bin/Raw\_FF\_et(mt).cc, FFcode/bin/OSSScorrection\_et(mt).cc, and FFcode/bin/Set1\_correction\_et(mt).cc (search for "Change here" and uncomment the WP you want)

To apply the FF:
----------------

Once the FF have been produced, they can be applied in your code as:


Include in the headers:

```
include "ComputeFF2018/FFcode/interface/ApplyFF.h"
```

Define all the functions:

```
   TFile frawff("ff_files/uncorrected_fakefactors_et.root");
   TF1* ff_qcd_0jet=(TF1*) frawff.Get("rawFF_et_qcd_0jet");
   TF1* ff_qcd_1jet=(TF1*) frawff.Get("rawFF_et_qcd_1jet");
   TF1* ff_w_0jet=(TF1*) frawff.Get("rawFF_et_w_0jet");
   TF1* ff_w_1jet=(TF1*) frawff.Get("rawFF_et_w_1jet");
   TF1* ff_tt_0jet=(TF1*) frawff.Get("mc_rawFF_et_tt");

   TFile fmvisclosure ("ff_files/FF_corrections_1.root");
   TF1* mvisclosure_qcd=(TF1*) fmvisclosure.Get("closure_mvis_et_qcd");
   TF1* mvisclosure_w=(TF1*) fmvisclosure.Get("closure_mvis_et_w");
   TF1* mvisclosure_tt=(TF1*) fmvisclosure.Get("closure_mvis_et_ttmc");

   TFile fosssclosure ("ff_files/FF_QCDcorrectionOSSS.root");
   TF1* osssclosure_qcd=(TF1*) fosssclosure.Get("closure_mvis_et_qcd");
   TF1* mtclosure_w=(TF1*) fosssclosure.Get("closure_mt_et_w");

```

Then compute the weight per event:


```
float my_fakefactor = get_ff(pt, mt, mvis, njets, frac_tt, frac_qcd, frac_w, ff_qcd_0jet, ff_qcd_1jet, ff_w_0jet, ff_w_1jet, ff_tt_0jet, mvisclosure_qcd, mvisclosure_w, mvisclosure_tt, mtclosure_w, osssclosure_qcd)
```

## Python Application Module
A python module is packaged with the repository, and should be available for use in python scripts applying fake factors. To use it, 
at the top of your python import the package with:

```import ComputeFF2018.FFcode.ApplyFF as ApplyFF``` 

To create the actual FF application tool, in the main body of your code do: 

```theFFApplicationTool = ApplyFF.FFApplicationTool("path/to/FFsForAYear/",channel)``` 

where channel is the channel you are creating FFs for, (either "et" or "mt" at the moment).
To retrieve fake factors, do : 

```event_fakefactor = theFFApplicationTool.get_ff(TauPt,TransverseMass,m_vis,njets,FracTT,FracQCD,FracW)```
