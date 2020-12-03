


count=`ls *.root | wc -l`
if [ $count == 89 ]
then 
    echo true
    rm mine_rootfile/*.root
    mv *.root mine_rootfile/
    echo "root files moved"

    cd mine_rootfile/
    hadd DY.root DYJetsToLL_*.root
    hadd DY1.root DY1JetsToLL_*.root
    hadd DY2.root DY2JetsToLL_*.root
    hadd DY3.root DY3JetsToLL_*.root
    hadd DY4.root DY4JetsToLL_*.root
    hadd TTTo2L2Nu.root TTTo2L2Nu_*.root
    hadd TTToHadronic.root TTToHadronic_*.root
    hadd TTToSemiLeptonic.root TTToSemiLeptonic_*.root
    hadd Data.root EGamma2018*.root
    cd ..
fi    
cd plotting_script/
bash zttPlots.sh
