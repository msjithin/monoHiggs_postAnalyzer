
declare -a rootFIles=("DY.root" "DY1.root" "DY2.root" "DY3.root" "DY4.root" "Data.root" "TTTo2L2Nu.root" "TTToHadronic.root" "TTToSemiLeptonic.root")
for FILE in "${rootFIles[@]}"
do
    if [ -a $FILE ]
    then
	echo "moving $FILE to theirs_rootfile/"
	mv $FILE theirs_rootfile/
    fi
done


count=`ls *.root | wc -l`
if [ $count == 70 ]
then 
    echo true
    rm mine_rootfile/*.root
    mv *.root mine_rootfile/
    echo "root files moved"

    cd mine_rootfile/
    hadd DY1.root DY1JetsToLL_M*.root
    hadd DY2.root DY2JetsToLL_M*.root
    hadd DY3.root DY3JetsToLL_M*.root
    hadd DY4.root DY4JetsToLL_M*.root
    hadd DY.root  DYJetsToLL_M*.root
    hadd Data.root SingleMuon*.root
    hadd TTTo2L2Nu.root TTTo2L2Nu_TuneCP5*.root
    hadd TTToHadronic.root TTToHadronic_TuneCP5_*.root
    hadd TTToSemiLeptonic.root TTToSemiLeptonic_TuneCP5_*.root


    cd ..
fi    
cd plotting_script/
bash zttPlots.sh
