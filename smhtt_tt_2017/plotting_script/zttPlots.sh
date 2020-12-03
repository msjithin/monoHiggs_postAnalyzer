set -e
# # keep track of the last executed command
# trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# # echo an error message before exiting
# trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

unset DISPLAY
echo "making plots ..........."

declare -a List_index=("_a" "_b" "_c" "_d" "_e" "_f" "_g" "_h" "_i")
declare -a List_names=("tau1Pt" "mT_eMet" "nJet" "mass1" "mass2" "genmatch1" "genmatch2" "met" "metPhi" "trigger" "tau2Pt" "tau1Eta" "tau2Eta" "tau1Phi" "tau2Phi" "deltaPhi" "deltaEta" "deltaR" "nEvents" "eventWeight" "higgsPt")

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
	python makeplot_ztt.py -bkg DY -name $hist -ch etau -xaxis $hist  -year 2017 &
	#python makeplot_ztt.py -bkg DY1 -name $hist -ch etau -xaxis $hist  -year 2017 &
        #python makeplot_ztt.py -bkg DY2 -name $hist -ch etau -xaxis $hist  -year 2017 &
	#python makeplot_ztt.py -bkg DY3 -name $hist -ch etau -xaxis $hist  -year 2017 &
	#python makeplot_ztt.py -bkg DY4 -name $hist -ch etau -xaxis $hist  -year 2017 &
	python makeplot_ztt.py -bkg TTTo2L2Nu -name $hist -ch etau -xaxis $hist  -year 2017 & 
	python makeplot_ztt.py -bkg TTToHadronic -name $hist -ch etau -xaxis $hist  -year 2017 &
	#python makeplot_ztt.py -bkg TTToSemiLeptonic -name $hist -ch etau -xaxis $hist  -year 2017 &
	#python makeplot_ztt.py -bkg DY_fbkg -name $hist -ch etau -xaxis $hist  -year 2017 &
	#wait
    done
    wait
done
wait
echo "Dy plots made"
for n in "${List_names[@]}"
do
    for i in "${List_index[@]}"
    do
	hist=$n$i
	#echo "$hist"
	python makeplot_ztt.py -bkg Data -name $hist -ch etau -xaxis $hist  -year 2017 &
	#python makeplot_ztt.py -bkg Data_fbkg -name $hist -ch etau -xaxis $hist  -year 2017
    done
    wait
done
wait
echo "Data plots made"

#sh combinne_to_pdf.sh

tar -czvf plots_2017.tar.gz plots_2017/
