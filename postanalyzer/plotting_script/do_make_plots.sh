set -e

declare -a plotList=("muPt" "tauPt" "higgsPt" "tot_TMass_full" "nJet" "visMass" "met" "metLongXaxis" "mT_muMet" "deltaR") 

declare -a indexList=("_"$1)


FILE=eventYield.csv
if [ -f "$FILE" ]; then
    echo "$FILE exists."
    rm eventYield.csv
fi
for i in "${indexList[@]}"
do
    for j in "${plotList[@]}"
    do 
	hist=$j$i
	python3 makeplot.py -name $hist -cat 0 -ch mutau -xaxis $hist -year 2017 --blindingRatio 1 
    done
    wait
done
