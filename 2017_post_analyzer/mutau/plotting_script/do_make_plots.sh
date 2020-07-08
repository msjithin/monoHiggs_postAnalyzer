

inFile=f_mutau_initial.root
#inFile=$1
declare -a plotList=("muPt_5" "muEta_5" "muPhi_5" "muDz_5" "muD0_5" "muonID_5" "relMuIso_5" "muCharge_5" "tauPt_5" "tauEta_5" "tauPhi_5" "tauIso_5" "tauDecayMode_5" "tauCharge_5" "tauAntiEle_5" "tauAntiMu_5" "deltaR_5" "higgsPt_5" "nJet_5" "visMass_5" "mT_muMet_5" "trigger_5" "muPt_6" "muEta_6" "muPhi_6" "muDz_6" "muD0_6" "muonID_6" "relMuIso_6" "muCharge_6" "tauPt_6" "tauEta_6" "tauPhi_6" "tauIso_6" "tauDecayMode_6" "tauCharge_6" "tauAntiEle_6" "tauAntiMu_6" "deltaR_6" "higgsPt_6" "nJet_6" "visMass_6" "mT_muMet_6" "trigger_6") 


for i in "${plotList[@]}"
do 
    echo "$i"
    #python makeplot_mutauh.py -in $inFile -name $i -cat 0 -ch mutau -xaxis $i
    python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name $i -cat 0 -ch mutau -xaxis $i -year 2017
done

#python makeplot_mutauh.py -in $inFile -name cutflow_n -cat 0 -ch mutau -xaxis "cutflow" -lY
python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name cutflow_n -cat 0 -ch mutau -xaxis "cutflow" -lY  -year 2017
