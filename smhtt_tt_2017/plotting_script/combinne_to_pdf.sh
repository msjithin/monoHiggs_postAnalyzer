
declare -a List_index=("j")
for i in "${List_index[@]}"
do
    inputName="plot_*_"$i"_etau_"
    echo "convert plots/Data/$inputName"Data.png" plots/Data_$i.pdf"
    convert plots/Data/$inputName"Data.png" plots/Data_$i.pdf &
    convert plots/DY/$inputName"DY.png" plots/DY_$i.pdf &
done

wait
echo "pdf files made!"
