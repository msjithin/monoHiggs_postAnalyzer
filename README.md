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
- All applied functions are in `src/object_functions.h`.

