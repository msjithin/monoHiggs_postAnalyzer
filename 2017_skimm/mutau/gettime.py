import os
directorylist = ['Out_DATA_03-05-2020_05-02', 'Out_DATA_06-05-2020_18-00', 'Out_DATA_06-05-2020_18-00', 'Out_DATA_29-04-2020_07-38', 'Out_MC_29-04-2020_07-37']


for directory in directorylist:
    for filename in os.listdir(directory):
        if filename.endswith(".out"):
            #print('In file {}'.format(filename))
            with open(directory+'/'+filename) as read_obj:
                for line in read_obj:
                    if 'Real time' in line:
                        line = line.rstrip()  # remove '\n' at end of line
                        print('{:<30s} {:>100s}'.format(filename, line))



