

inFile=f_mutau_initial.root
#inFile=$1
declare -a plotList=("muPt_Plus_5" "muEta_Plus_5" "muPhi_Plus_5" "muDz_Plus_5" "muD0_Plus_5" "muonID_Plus_5" "relMuIso_Plus_5" "muCharge_Plus_5" "muPt_Minus_5" "muEta_Minus_5" "muPhi_Minus_5" "muDz_Minus_5" "muD0_Minus_5" "muonID_Minus_5" "relMuIso_Minus_5" "muCharge_Minus_5" "muMuMass_5" ) 


for i in "${plotList[@]}"
do 
    echo "$i"
    j=${i::-2}
    #echo "    python makeplot_mutauh.py -in $inFile -name $i -cat 0 -ch mutau -xaxis $i "
    python makeplot_mutauh.py -in $inFile -name $i -cat 0 -ch mutau -xaxis $j
    python makeplot_mutauh.py -in $inFile -name $i -cat 0 -ch mutau -xaxis $j -lY

done


#for i in "${plotList_1[@]}"
#do
#    echo "$i"
#    python makeplot_mutauh.py -in $inFile -name $i -cat 0 -ch mutau -xaxis $i
#done

#python makeplot_mutauh.py -in $inFile -name cutflow_n -cat 0 -ch mutau -xaxis "cutflow" -lY
#python makeplot_mutauh.py -in $inFile -name cutflow_Htt -cat 0 -ch mutau -xaxis "cutflow" -lY

#python makeplot_mutauh.py -in $inFile -lY -name Higgs_pt -cat 0 -ch mutau
