

myEvents=[]
fp_mine=  open('mine_rootfile/dyjets.txt', 'r')
for x in fp_mine:
    split_string = x.split()
    #print(split_string[0])
    myEvents.append(split_string[0])

their_events=[]
fp_theirs= open('theirs_rootfile/eventAnalysis.txt', 'r')
for x in fp_theirs:
    split_string = x.split()
    #print(split_string[0])
    their_events.append(split_string[0])

file2 = open("output.txt","w+") 
rogue_event_list=[]
good_event_list=[]
count=1
for i in myEvents:
    if i not in their_events:
        rogue_event_list.append(i)
        #file2.write(i+'\n')
    if i in their_events:
        print(count, i)
        good_event_list.append(i)
        count+=1
file2.write("rogue_event_list = "+str(rogue_event_list))
div_str="* \n "
file2.write(div_str*10)
file2.write("good_event_list = "+str(good_event_list))
