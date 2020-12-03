set -e
# # keep track of the last executed command
# trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# # echo an error message before exiting
# trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

unset DISPLAY
echo "making plots ..........."

declare -a List_index=("_a" "_b" "_c" "_d" "_e" "_f" "_g" "_h" "_i" "_j")
declare -a List_names=("muPt" "mT_muMet" "tauDecayMode" "nJet" "mass1" "mass2" "genmatch1" "genmatch2" "met" "metPhi" "trigger" "tauPt" "muEta" "tauEta" "muPhi" "tauPhi" "deltaPhi" "deltaEta" "deltaR" "nEvents" "eventWeight" "higgsPt" "fakefactor")

if [ -f "eventYield.csv" ]; then
    rm eventYield.csv
fi
#rm plots/DY/*.png
#rm plots/Data/*.png
for n in "${List_names[@]}"
do
    for i in "${List_index[@]}"
    do
	hist=$n$i
	#echo "$hist"
	python makeplot_ztt.py -bkg DY -name $hist -ch mutau -xaxis $hist  -year 2017 & 
	python makeplot_ztt.py -bkg DY1 -name $hist -ch mutau -xaxis $hist  -year 2017 & 
        python makeplot_ztt.py -bkg DY2 -name $hist -ch mutau -xaxis $hist  -year 2017 & 
	python makeplot_ztt.py -bkg DY3 -name $hist -ch mutau -xaxis $hist  -year 2017 & 
	python makeplot_ztt.py -bkg DY4 -name $hist -ch mutau -xaxis $hist  -year 2017 & 
	python makeplot_ztt.py -bkg TTTo2L2Nu -name $hist -ch mutau -xaxis $hist  -year 2017 & 
	python makeplot_ztt.py -bkg TTToHadronic -name $hist -ch mutau -xaxis $hist  -year 2017 & 
	python makeplot_ztt.py -bkg TTToSemiLeptonic -name $hist -ch mutau -xaxis $hist  -year 2017 & 
	#python makeplot_ztt.py -bkg DY_fbkg -name $hist -ch mutau -xaxis $hist  -year 2017 & 
	wait
    done
    #wait
done
wait
echo "Dy plots made"
for n in "${List_names[@]}"
do
    for i in "${List_index[@]}"
    do
	hist=$n$i
	#echo "$hist"
	python makeplot_ztt.py -bkg Data -name $hist -ch mutau -xaxis $hist -year 2017 & 
	#python makeplot_ztt.py -bkg Data_fbkg -name $hist -ch mutau -xaxis $hist  -year 2017  &
    done
    #wait
done
wait
echo "Data plots made"

echo "All processes done!"
sh combinne_to_pdf.sh

tar -czvf plots_2017.tar.gz plots_2017/
