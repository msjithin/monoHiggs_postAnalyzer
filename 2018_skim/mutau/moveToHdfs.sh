

#rsync -zvh --progress --remove-source-files DY1JetsToLL_00.root /hdfs/store/user/jmadhusu/2018_skimmed/

#rsync -zvh --progress --remove-source-files *.root /hdfs/store/user/jmadhusu/2018_skimmed/
echo -n """
 moving files to /hdfs/store/user/jmadhusu/2018_skimmed/mutau/
"""

for jobFile in `ls output`; do  
 echo -n """  
 ******************************
 moving  ${jobFile}
"""   
 mv output/${jobFile} /hdfs/store/user/jmadhusu/2018_skimmed/mutau/
 

done   
