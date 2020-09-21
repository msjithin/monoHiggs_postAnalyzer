

count=`ls *.root | wc -l`
if [ $count == 52 ]
then 
    echo true
    rm mine_rootfile/*.root
    mv *.root mine_rootfile/
    echo "root files moved"

    cd mine_rootfile/
    hadd DY.root DYJetsToLL_*.root
    hadd Data.root EGamma2018*.root
    cd ..
fi
cd plotting_script/
bash zttPlots.sh
