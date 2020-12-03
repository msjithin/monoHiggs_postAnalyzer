
declare -a List_index=("j")
for i in "${List_index[@]}"
do
    inputName="plot_*_"$i"_etau_"
    echo "convert plots_2017/Data/$inputName"Data.png" plots_2017/Data_$i.pdf"
    convert plots_2017/Data/$inputName"Data.png" plots_2017/Data_$i.pdf &
    convert plots_2017/DY/$inputName"DY.png" plots_2017/DY_$i.pdf &
done

wait
echo "pdf files made!"
