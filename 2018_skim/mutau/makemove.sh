

for (( i=0; i<20; i++)); do

    date
    echo -n """
    ${i}
    moving files to /hdfs/store/user/jmadhusu/2018_skimmed/mutau/
    """
    mv *.root output/
    for jobFile in `ls output`; do  
	echo -n """  
        ******************************
         moving  ${jobFile}
        """   
	mv output/${jobFile} /hdfs/store/user/jmadhusu/2018_skimmed/mutau/
	
    done   
    
    
    echo "Sleep for 30 min"
    
    sleep 1800
    
done

exit 0;
