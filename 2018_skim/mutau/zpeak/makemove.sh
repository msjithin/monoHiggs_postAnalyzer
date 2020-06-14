
#rsync -zvh --progress --remove-source-files DY1JetsToLL_00.root /hdfs/store/user/jmadhusu/2018_skimmed/

#rsync -zvh --progress --remove-source-files *.root /hdfs/store/user/jmadhusu/2018_skimmed/
for (( i=0; i<24; i++)); do

    date
    echo -n """
    ${i}
    moving files to /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/
    """
    mv *.root output/
    for jobFile in `ls output`; do  
	echo -n """  
        ******************************
         moving  ${jobFile}
        """   
	mv output/${jobFile} /hdfs/store/user/jmadhusu/2018_skimmed/muSelections/ 
	
    done   
    
    
    echo "Sleep for 15 min"
    
    sleep 900
    
done

exit 0;
