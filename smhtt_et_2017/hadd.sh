


count=`ls *.root | wc -l`
if [ $count == 33 ]
then 
    echo true
    rm mine_rootfile/*.root
    mv *.root mine_rootfile/
    echo "root files moved"

    cd mine_rootfile/
    hadd DY.root DYJetsToLL_M-50*.root
    hadd Data.root SingleElectron_*.root
    cd ..
fi    
cd plotting_script/
bash zttPlots.sh
