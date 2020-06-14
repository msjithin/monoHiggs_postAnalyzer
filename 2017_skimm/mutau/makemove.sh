

for (( i=0; i<10; i++)); do

    date
    echo -n """
    ${i}
    moving files to /hdfs/store/user/jmadhusu/2017_skimmed/mutau/
    """
    mv *.root output/
    for jobFile in `ls output`; do  
	echo -n """  
        ******************************
         moving  ${jobFile}
        """   
	mv output/${jobFile} /hdfs/store/user/jmadhusu/2017_skimmed/mutau/
	
    done   
    
    
    echo "Sleep for 30 min"
    
    sleep 1800
    
done

exit 0;
