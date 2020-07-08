

for (( i=0; i<40; i++)); do

    date
    echo -n """
    ${i}
    moving files to /hdfs/store/user/jmadhusu/2017_skimmed/etau
    """
    mv *.root output/
    for jobFile in `ls output`; do  
	echo -n """  
        ******************************
         moving  ${jobFile}
        """   
	mv output/${jobFile} /hdfs/store/user/jmadhusu/2017_skimmed/etau/
	
    done   
    
    
    echo "Sleep for 15 min"
    
    sleep 900
    
done

exit 0;
