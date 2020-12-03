

inFile=f_mutau_initial.root

declare -a plotList_a=("elePt_a" "elePt_b" "elePt_c" "elePt_d" "elePt_e" "elePt_f" "elePt_g" "elePt_h" "elePt_i" "elePt_j")

for i in "${plotList_a[@]}"
do
   echo "$i"
   #python ~/monoHiggs_2018_wDnn/CMSSW_10_2_18/src/analysis/MacrosAndScripts/makeplot.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i -year 2017
   python makeplot_ztt.py -in $inFile -name $i -cat 0 -ch etau -xaxis $i -year 2017
done

