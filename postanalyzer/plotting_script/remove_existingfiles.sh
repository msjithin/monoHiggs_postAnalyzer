


agg_file(){
    var_name=$1
    if [[ $# -eq 0 ]] ; then
        echo 'no arguments given, which files to remove?'
        exit 0
    fi
    if [ "$(ls -A *${var_name}*.root)" ]; then                                                                                                   
        echo "Take action root files exits"                                                                                          
        rm *${var_name}*.root                                                                                                                 
    fi                                                                                                                               
}


scaled_files(){
    var_name=$1
    if [[ $# -eq 0 ]] ; then
        echo 'no arguments given, which files to remove?'
        exit 0
    fi
    outDIR="sample"                                                                                                                  
    if [ -d "$outDIR" ]; then                                                                                                        
        echo "$outDIR exists"                                                                                                        
        if [ "$(ls -A $outDIR)" ]; then                                                                                              
            echo "Take action $outDIR is not Empty .... removing existing files ....."  
            if [ "$(ls -A *${var_name}*.root)" ]; then
            echo "removing $outDIR/*${var_name}*.root  "
                rm $outDIR/*${var_name}*.root   
            fi                                                                                                                                             
        else                                                                                                                         
            echo " $outDIR is Empty"                                                                                                 
        fi                                                                                                                           
    else                                                                                                                             
        echo "$outDIR created"                                                                                                       
        mkdir $outDIR                                                                                                                
    fi                  
}


"$@"