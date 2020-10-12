


myFile = open('../out_sf', 'r')
myLines = myFile.readlines()
outstr="{"
count=0
for line in myLines:
    myList=line.strip().split(",")
    #print myList
    if  myList[0]=='event':
        #print myList
        notFoundEvent=True
        with open('../out_sf_smhtt') as f:
            for thline in f:
                thList=thline.strip().split(",")
                if thList[0]=='event':
                    #print(myList[1], thList[1])
                    if( myList[1]==thList[1]):
                        notFoundEvent=False
        if notFoundEvent:
            count=count+1
            outstr=outstr+myList[1]+","

outstr=outstr+"}"

print(outstr)
print(count)
