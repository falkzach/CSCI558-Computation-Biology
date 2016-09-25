#!/usr/local/bin/python3

myList = [8,1,5,1,0,4,5,10,1]
myList.sort()
newList = []
last = None
for x in myList:
    if x is not last:
        newList.append(x)
    last = x

print(newList)
