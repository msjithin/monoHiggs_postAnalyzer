# monoHiggs to tautau analysis code

# mutau channel

All scripts are in directory `src`.

I use `runAnalyzer.sh` for running interactively and locally.
Use `python step1_makeSubmitFile.py` to genertae submit file for all MC and data.  
This will generate `submit_condor.sh`

To submit all jobs:  
`bash submit_condor.sh`



The main scripts are `src/mutau_analyzer.h and .C`.
Selections are applied in `src/selections.h`

Workflow:  
- First select muons, taus, boosted taus.
- If atleast one muon and one tau is found, use `get_index` function to get index of muon and tau/boostedtau.
- mu-tau dr > 0.5 ==> deeptau isolations are applied.  
mu-tau dr < 0.5 ==> boostedTau isolations are applied.
- `setMyEleTau()` function set the 4-vector of muon, tau, and related variables to be used for further selections. This function is in mutau_analyzer.h .
- All applied functions are in `src/object_functions.h`

There is a processing step once all the root files are generated. 
- Move all files to `output` directory.
- execute  `bash execute_all.sh ` 
- execute ` bash do_make_plots.sh ` to generate plots.

[The plots can be viewd here : display_plots.html](https://htmlpreview.github.io/?https://github.com/msjithin/monoHiggs_postAnalyzer/blob/mutau_2017/display_plots.html)


Selections:
- Trigger : HLT_PFMET120_PFMHT120_IDTight_PFHT60_v, HLT_PFMET120_PFMHT120_IDTight_v
- Muons : pt > 20 , eta < 2.4
    - medium muon id
- Taus :  pt> 30 , eta < 2.3 
    - boosted Medium Isolation
    - boosted Loose electron rejection
    - boosted loose muon rejection
- mu-tau opposite charge
- mu-tau deltaR < 0.5
- Higgs pt > 65
- mu-tau visible mass < 125
- MET > 105