set -e

inFile=f_tautau_initial.root
#inFile=$1
declare -a plotList=("tau1Pt" "tau1Eta" "tau1Phi" "tau2Pt" "tau2Eta" "tau2Phi" "tau1Charge" "tau1DecayMode" "tau2DecayMode" "tau2Charge" "tau1Iso" "tau2Iso" "tau1Charge" "tau2Charge" "deltaR" "deltaPhi" "deltaEta" "nJet" "met" "metLongXaxis" "metPhi" "mT_tauMet" "visMass" "higgsPt" "trigger") 
declare -a indexList=("_3" "_4" "_5")

if [ -f "eventYield.csv" ]; then
    echo "removing eventYield file"
    rm eventYield.csv
fi


for i in "${indexList[@]}"
do
    for j in "${plotList[@]}"
    do 
	hist=$j$i
	#echo "$hist"
	#python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i -year 2017 
	python makeplot.py -in $inFile -name $hist -cat 0 -ch tautau -xaxis $hist -year 2017 &
	python makeplot_fbkg.py -in $inFile -name $hist -cat 0 -ch tautau -xaxis $hist -year 2017 &
    done
    wait
done

wait
python makeplot.py -in $inFile -name cutflow_n -cat 0 -ch tautau -xaxis "cutflow" -lY  -year 2017 
