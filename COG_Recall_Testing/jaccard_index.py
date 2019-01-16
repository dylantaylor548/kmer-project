import os

def text2list(file_path):
    textlist = []

    with open(file_path) as opened_file:
        for line in opened_file:
            line = line.strip()
            textlist.append(line)
    return textlist


# Takes two sets as inputs
def jaccard(a, b):
    c = a.intersection(b)
    print("There are " + str(len(c)) + " kmers shared between the two lists.")
    return float(len(c)) / (len(a) + len(b) - len(c))


from itertools import combinations

# change if you choose a different number of partitions

directory = input("Input the path to the directory containing your representative list files: ")

listofreplists = {}
i=1

for file_name in os.listdir(directory):
    if file_name.endswith('.txt'):
        file_path = (directory + '/' + file_name)
        replist = text2list(file_path)
        print(len(replist))
        listofreplists[str(i)] = replist
        i += 1

partitions= []
for x in range(1,i):
    partitions.append(x)

comb = combinations(partitions,2)
for (x,y) in list(comb):
    xset = set(listofreplists[str(x)])
    yset = set(listofreplists[str(y)])
    temp = jaccard(xset,yset)
    print("Jaccard Index for lists: " +str(x) + ", " + str(y) + " is " +str(temp))
