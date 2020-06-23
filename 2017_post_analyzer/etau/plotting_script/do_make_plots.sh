

inFile=f_mutau_initial.root
#inFile=$1
declare -a plotList=("elePt_5" "eleEta_5" "elePhi_5" "eleDz_5" "eleD0_5" "electronID_5" "relEleIso_5" "eleCharge_5" "tauPt_5" "tauEta_5" "tauPhi_5" "tauIso_5" "tauDecayMode_5" "tauCharge_5" "tauAntiEle_5" "tauAntiMu_5" "deltaR_5" "higgsPt_5" "nJet_5" "visMass_5" "mT_muMet_5" "elePt_6" "eleEta_6" "elePhi_6" "eleDz_6" "eleD0_6" "electronID_6" "relEleIso_6" "eleCharge_6" "tauPt_6" "tauEta_6" "tauPhi_6" "tauIso_6" "tauDecayMode_6" "tauCharge_6" "tauAntiEle_6" "tauAntiMu_6" "deltaR_6" "higgsPt_6" "nJet_6" "visMass_6" "mT_muMet_6") 


for i in "${plotList[@]}"
do 
    echo "$i"
    python makeplot_mutauh.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i
done


#for i in "${plotList_1[@]}"
#do
#    echo "$i"
#    python makeplot_mutauh.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i
#done

python makeplot_mutauh.py -in $inFile -name cutflow_n -cat 0 -ch etau -xaxis "cutflow" -lY
#python makeplot_mutauh.py -in $inFile -name cutflow_Htt -cat 0 -ch etau -xaxis "cutflow" -lY
#python makeplot_mutauh.py -in $inFile -lY -name Higgs_pt -cat 0 -ch mutau
