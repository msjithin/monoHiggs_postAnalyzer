import os
directorylist = ['Out_DATA_29-04-2020_10-51', 'Out_DATA_30-04-2020_00-51', 'Out_MC_29-04-2020_10-52']


for directory in directorylist:
    for filename in os.listdir(directory):
        if filename.endswith(".out"):
            #print('In file {}'.format(filename))
            with open(directory+'/'+filename) as read_obj:
                for line in read_obj:
                    if 'Real time' in line:
                        line = line.rstrip()  # remove '\n' at end of line
                        print('{:<30s} {:>100s}'.format(filename, line))



