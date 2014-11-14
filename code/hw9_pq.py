# Written in Python 3.3.2
# Author:  Rudy Laprade

import queue, time

def main():
    prev_t = time.clock()
    for i in range(1, 17):
        input_file = "../Dataset/reads{}.txt".format(i)
        dna = reconstruct(input_file) + '\n'
        with open("../output{}.txt".format(i), "w") as output:
            output.write(dna)
            if i < 6:
                with open("../Dataset/answer{}.txt".format(i), "r") as answer:
                    if dna == answer.read():
                        print("Output correct for #{}".format(i))
                    else:
                        print("Incorrect output for #{}".format(i))
        curr_t = time.clock()
        print("Time to complete {}:  {} seconds".format(i, curr_t-prev_t))
        prev_t = curr_t
    print("Total running time: {} min, {} sec".format(curr_t//60, curr_t%60))

def reconstruct(file):
    pq = queue.PriorityQueue()
    with open(file, "r") as f:
        strand = {}
        key = 0
        for line in f:
            strand[key] = line[:(len(line)-1)]
            key += 1
        
    # Calculate overlaps and put pairs in a priority queue prioritizing largest overlap
    for i in range(len(strand)):
        for j in range(len(strand)):
            if i != j:
                pq.put((-overlap(strand[i],strand[j]), (i,j)))
                
    print("Time to compute initial overlaps {}".format(time.clock()))
                
    while not pq.empty():
        over, pair = pq.get(False)
        i, j = pair
        over = abs(over)
        # If both strands valid, merge them and add merged strand to end of list
        if strand[i] is not None and strand[j] is not None:
        # if i in strand and j in strand:
            newStr = combine(strand[i], strand[j], over)
            strand[i], strand[j] = None, None # invalidate previous strands
            strand[key] = newStr
            k = key
            key += 1
            # Calculate overlaps between merged strand and other strands
            for index in range(k):
                if strand[index] != None:
                    pq.put((-overlap(strand[k],strand[index]), (k,index)))
                    pq.put((-overlap(strand[index],strand[k]), (index,k)))
    return strand[k]
    
def overlap(str1, str2):
    for i in range(len(str1)):
        if len(str1)-i > len(str2):
            if str1[i:i+len(str2)] == str2:
                return len(str2)
        elif str1[i:] == str2[:len(str1)-i]:
            return len(str1)-i
    return 0
    
def combine(str1, str2, overlap):
    return str1 + str2[overlap:]

if __name__ == "__main__":
    main()