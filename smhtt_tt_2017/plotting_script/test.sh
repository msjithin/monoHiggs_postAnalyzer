
unset DISPLAY

echo "making plots ..........."
declare -a List_index=("_a" "_b" "_c" "_d" "_e" "_f" "_g" "_h" "_i" "_j")
declare -a List_names=("met" "metPhi" "trigger")


for n in "${List_names[@]}"
do
    for i in "${List_index[@]}"
    do
	hist=$n$i
	echo "$hist"

    done
done

wait
echo "All processes done!"
