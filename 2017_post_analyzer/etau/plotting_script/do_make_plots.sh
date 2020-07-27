

inFile=f_etau_initial.root
#inFile=$1
declare -a plotList=("elePt_5" "eleEta_5" "elePhi_5" "eleDz_5" "eleD0_5" "electronID_5" "relEleIso_5" "eleCharge_5" "tauPt_5" "tauEta_5" "tauPhi_5" "tauIso_5" "tauDecayMode_5" "tauCharge_5" "tauAntiEle_5" "tauAntiMu_5" "deltaR_5" "higgsPt_5" "nJet_5" "visMass_5" "mT_eleMet_5" "trigger_5" "genMatch_5" "elePt_6" "eleEta_6" "elePhi_6" "eleDz_6" "eleD0_6" "electronID_6" "relEleIso_6" "eleCharge_6" "tauPt_6" "tauEta_6" "tauPhi_6" "tauIso_6" "tauDecayMode_6" "tauCharge_6" "tauAntiEle_6" "tauAntiMu_6" "deltaR_6" "higgsPt_6" "nJet_6" "visMass_6" "mT_eleMet_6" "trigger_6" "genMatch_6" "met_5" "met_6") 

declare -a plotList_a=("elePt_a" "elePt_b" "elePt_c" "elePt_d" "elePt_e" "elePt_f" "elePt_g" "elePt_h" "elePt_i" "elePt_j")

for i in "${plotList[@]}"
do 
    echo "$i"
    #python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i -year 2017 
    python makeplot.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i -year 2017
done


for i in "${plotList_a[@]}"
do
   echo "$i"
   #python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i -year 2017
   python makeplot.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i -year 2017
   
done
for i in "${plotList_a[@]}"
do
   echo "$i"
   #python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i -year 2017
   python makeplot_ztt.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i -year 2017
done

#python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name cutflow_n -cat 0 -ch etau -xaxis "cutflow" -lY  -year 2017
python makeplot.py -in $inFile -name cutflow_n -cat 0 -ch etau -xaxis "cutflow" -lY  -year 2017
