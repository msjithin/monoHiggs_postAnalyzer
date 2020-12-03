


myFile = open('out_mono', 'r')
myLines = myFile.readlines()
outstr="{"
count=0
for line in myLines:
    myList=line.strip().split("\t")
    #print myList
    if  myList[0]=='event':
        #print myList
        notFoundEvent=True
        with open('out_smhtt') as f:
            for thline in f:
                thList=thline.strip().split("\t")
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
