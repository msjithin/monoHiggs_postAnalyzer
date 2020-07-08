# Install
This has been tested on CMSSW_7_6_3  

`git clone https://github.com/CMS-HTT/Jet2TauFakes.git HTTutilities/Jet2TauFakes`  
`scram b -j4`   

# Tests
`cd HTTutilities/Jet2TauFakes/test`   

## C++
`root`   
 `.x loadLibrary.C`   
 `.L test2.C`   
 `test2()`   

## Python
Test direct access to fake-factor objects:  
`python test.py`  
Test python utilities for more user-friendly access:  
`python test_utilities.py`  
Test systematic definitions:  
`python test_sys.py` 


