set -e



# python main.py -hist 4,7,19,20,21,22,51,26,30 -idx 4 -ch mutau -y 2017
# bash do_make_plots.sh 4

# python main.py -hist 4,7,19,20,21,22,51,26,30 -idx 5 -ch mutau -y 2017
# bash do_make_plots.sh 5

# # python main.py -hist 4,7,19,20,21,22,51,26,30 -idx 6 -ch mutau -y 2017
# # bash do_make_plots.sh 6

# python main.py -hist 4,7,19,20,21,22,51,26,30 -idx 7 -ch mutau -y 2017
# bash do_make_plots.sh 7

# python main.py -hist 4,7,19,20,21,22,51,26,30 -idx 8 -ch mutau -y 2017
# bash do_make_plots.sh 8

python main.py -hist 4,7,19,20,21,22,51,26,30,32 -idx 9 -ch mutau -y 2017
bash do_make_plots.sh 9

# for j in 5 6 7 8 9
#     do 
#     for i in 4 7 19 20 21 22 51 26 30
#     do 
#         python main.py -hist $i -idx $j -ch mutau -y 2017 &
#     done
#     wait
#     bash do_make_plots.sh $j
# done

