#!/bin/bash
set -e  # exit when any command fails


./rootcom smhet_2018 executable_smhtt_et
./executable_smhtt_et theirs_rootfile/fromCecile_DY.root theirs_rootfile/DY.root -1 1000 2018 MC DY_v1 
./executable_smhtt_et theirs_rootfile/fromCecile_Data.root theirs_rootfile/Data.root -1 1000 2018 DATA Data 
