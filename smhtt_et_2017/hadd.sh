


count=`ls *.root | wc -l`
if [ $count == 66 ]
then 
    echo true
    rm mine_rootfile/*.root
    mv *.root mine_rootfile/
    echo "root files moved"

    cd mine_rootfile/
    hadd DY.root DYJetsToLL_M-50*00.root DYJetsToLL_M-50*01.root DYJetsToLL_M-50*02.root DYJetsToLL_M-50*03.root DYJetsToLL_M-50*04.root DYJetsToLL_M-50*05.root
    hadd DY_fbkg.root DYJetsToLL_M-50*_fbkg_.root
    hadd Data.root SingleElectron_*00.root  SingleElectron_*01.root SingleElectron_*02.root SingleElectron_*03.root SingleElectron_*04.root SingleElectron_*05.root
    hadd Data_fbkg.root SingleElectron_*_fbkg_.root
    cd ..
fi    
cd plotting_script/
bash zttPlots.sh
