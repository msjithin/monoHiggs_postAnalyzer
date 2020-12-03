set -e

inFile=f_mutau_initial.root
#inFile=$1
declare -a plotList=("muPt" "muEta" "muPhi" "muDz" "muD0" "muonID" "relMuIso" "muCharge" "tauPt" "tauEta" "tauPhi" "tauIso" "tauDecayMode" "tauCharge" "tauAntiEle" "tauAntiMu" "deltaR" "higgsPt" "nJet" "visMass" "mT_muMet" "trigger" "genMatch" "met" "metPhi" "deltaPhi" "deltaEta" "metLongXaxis") 
declare -a indexList=("_5" "_6")

if [ -f "eventYield.csv" ]; then
    rm eventYield.csv
fi
for i in "${indexList[@]}"
do
    for j in "${plotList[@]}"
    do
        hist=$j$i
        python makeplot.py -in $inFile -name $hist -cat 0 -ch mutau -xaxis $hist -year 2017 &
        #python makeplot_fbkg.py -in $inFile -name $hist -cat 0 -ch etau -xaxis $hist -year 2017 &
    done
    wait
done

wait
python makeplot.py -in $inFile -name cutflow_n -cat 0 -ch mutau -xaxis "cutflow" -lY  -year 2017




# for i in "${plotList[@]}"
# do 
#     echo "$i"
#     #python makeplot_mutauh.py -in $inFile -name $i -cat 0 -ch mutau -xaxis $i
#     python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name $i -cat 0 -ch mutau -xaxis $i -year 2017
# done

# #python makeplot_mutauh.py -in $inFile -name cutflow_n -cat 0 -ch mutau -xaxis "cutflow" -lY
# python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name cutflow_n -cat 0 -ch mutau -xaxis "cutflow" -lY  -year 2017
