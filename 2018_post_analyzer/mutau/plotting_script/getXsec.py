import os

filelist=os.listdir("xslogs")

for fileIn in filelist :
    sampleName=fileIn.strip('xsec_').strip('.log')
    print(sampleName)
    with open('xslogs/'+fileIn, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            if "final cross section" in line:
                print(line)

