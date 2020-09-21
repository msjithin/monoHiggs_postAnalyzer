
unset DISPLAY
echo "making plots ..........."

declare -a List_index=("_a" "_b" "_c" "_d" "_e" "_f" "_g" "_h" "_i" "_j")
declare -a List_names=("elePt" "mT_eMet" "tauDecayMode" "nJet" "decayModeFinding" "mass1" "mass2" "genmatch1" "genmatch2" "met" "metPhi" "trigger" "tauPt" "eleEta" "tauEta" "elePhi" "tauPhi" "deltaPhi" "deltaEta" "deltaR")


for n in "${List_names[@]}"
do
    for i in "${List_index[@]}"
    do
	hist=$n$i
	#echo "$hist"
	python makeplot_ztt.py -bkg DY -name $hist -ch etau -xaxis $hist  -year 2018 &
	python makeplot_ztt.py -bkg Data -name $hist -ch etau -xaxis $hist  -year 2018 & 
    done
done

wait
echo "All processes done!"

sh combinne_to_pdf.sh
