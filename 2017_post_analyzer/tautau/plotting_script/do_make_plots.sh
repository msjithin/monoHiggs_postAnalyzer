

inFile=f_mutau_initial.root
#inFile=$1
declare -a plotList=("tau1Pt_5" "tau1Eta_5" "tau1Phi_5" "tau1Iso_5" "tau1DecayMode_5" "tau1Charge_5" "tau1AntiEle_5" "tau1AntiMu_5" "tau2Pt_5" "tau2Eta_5" "tau2Phi_5" "tau2Iso_5" "tau2DecayMode_5" "tau2Charge_5" "tau2AntiEle_5" "tau2AntiMu_5" "deltaR_5" "higgsPt_5" "nJet_5" "visMass_5" "mT_muMet_5" "trigger_5" "tau1Pt_6" "tau1Eta_6" "tau1Phi_6" "tau1Iso_6" "tau1DecayMode_6" "tau1Charge_6" "tau1AntiEle_6" "tau1AntiMu_6" "tau2Pt_6" "tau2Eta_6" "tau2Phi_6" "tau2Iso_6" "tau2DecayMode_6" "tau2Charge_6" "tau2AntiEle_6" "tau2AntiMu_6" "deltaR_6" "higgsPt_6" "nJet_6" "visMass_6" "mT_muMet_6" "trigger_6") 


for i in "${plotList[@]}"
do 
    echo "$i"
    python makeplot.py -in $inFile -name $i -cat 0 -ch tautau -xaxis $i -year 2017
    #python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name $i -cat 0 -ch mutau -xaxis $i -year 2017
done

python makeplot.py -in $inFile -name cutflow_n -cat 0 -ch mutau -xaxis "cutflow" -lY -year 2017
#python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name cutflow_n -cat 0 -ch mutau -xaxis "cutflow" -lY  -year 2017
