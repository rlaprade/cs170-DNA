# Written in Python 3.3.2
# Author:  Rudy Laprade

import heapq, time

kmp_time = 0

def main():
    prev_t = time.clock()
    for i in range(1, 17):
        input_file = "./../Dataset/reads{}.txt".format(i)
        dna = reconstruct(input_file) + '\n'
        with open("./../output{}.txt".format(i), "w") as output:
            output.write(dna)
            if i < 6:
                with open("./../Dataset/answer{}.txt".format(i), "r") as answer:
                    if dna == answer.read():
                        print("Output correct for #{}".format(i))
                    else:
                        print("Incorrect output for #{}".format(i))
        curr_t = time.clock()
        print("Time to complete {}:  {} seconds \n".format(i, curr_t-prev_t))
        prev_t = curr_t
    print("kmp time {}".format(kmp_time))
    print("Total running time: {} min, {} sec".format(curr_t//60, curr_t%60))

def reconstruct(file):
    dfa_dict = {}       
    over_dict = {}
    def overlap(str1, str2):
        if str2 in dfa_dict:
            dfa = dfa_dict[str2]
        else:
            dfa = kmp(str2)
            dfa_dict[str2] = dfa
        len2 = len(str2)
        j = 0
        for c in str1:
            if j == len2:
                return j
            j = dfa[j][c]
        over_dict[(str1, str2)] = j
        return j    

    pq = []
    with open(file, "r") as f:
        strand = {}
        key = 0
        for line in f:
            strand[key] = line[:(len(line)-1)]
            key += 1
        
    # Calculate overlaps and put pairs in a priority queue prioritizing largest overlap
    begin = time.clock()
    for i in range(len(strand)):
        for j in range(len(strand)):
            if i != j and strand[i] != None and strand[j] != None and overlap(strand[i],strand[j]) == len(strand[j]):
                strand[j] = None
    # best_leading = {}
    # best_ending = {}
    # for i in range(len(strand)):
        # if strand[i] != None:
            # for j in range(len(strand)):
                # if strand[j] != None and i != j:
                    # if i not in best_leading or i not in best_ending:
                        # best_leading[i] = j
                        # best_ending[i] = j
                    # else:
                        # if overlap(strand[j],strand[i]) > overlap(strand[best_leading[i]], strand[i]):
                            # best_leading[i] = j
                        # if overlap(strand[i],strand[j]) > overlap(strand[i], strand[best_ending[i]]):
                            # best_ending[i] = j
    # for i in range(len(strand)):
        # if strand[i] != None:
            # heapq.heappush(pq, (-overlap(strand[best_leading[i]],strand[i]), (best_leading[i],i)))
            # heapq.heappush(pq, (-overlap(strand[i],strand[best_ending[i]]), (i,best_ending[i])))
    for i in range(len(strand)):
        if strand[i] != None:
            for j in range(len(strand)):
                if strand[j] != None and i != j:
                    heapq.heappush(pq, (-overlap(strand[i],strand[j]), (i,j)))
    print("Time to compute initial overlaps {}".format(time.clock()-begin))
                
    begin = time.clock()          
    while pq:
        over, pair = heapq.heappop(pq)
        i, j = pair
        over = abs(over)
        # If both strands valid, merge them and add merged strand to end of strand data structure
        if strand[i] is not None and strand[j] is not None:
            if over == len(strand[j]):
                strand[j] = None
                continue
            newStr = combine(strand[i], strand[j], over)
            strand[key] = newStr
            # k_len = len(newStr)
            l_str = strand[i]
            r_str = strand[j]
            # Calculate overlaps between merged strand and other strands
            for index in range(key):
                other = strand[index]
                if other != None:
                    rem = other[overlap(l_str,other):]
                    if rem == r_str[over:over+len(rem)]:
                        heapq.heappush(pq, (-len(other),(key,index)))
                    else:
                        heapq.heappush(pq, (-overlap(other,l_str),(index,key)))
                        heapq.heappush(pq, (-overlap(r_str, other),(key,index)))
            key += 1
            strand[i], strand[j] = None, None # invalidate previous strands
    print("Time to compute other part {}".format(time.clock()-begin))
    return strand[key-1]
    

def kmp(pat):
    global kmp_time
    start_time = time.clock()
    alpha = ('A', 'C', 'G', 'T')
    k = len(pat)
    dfa = [{} for _ in range(k)]
    for char in alpha: dfa[0][char] = 0
    dfa[0][pat[0]] = 1
    X = 0
    for j in range(1, k):
        for char in alpha:
            dfa[j][char] = dfa[X][char]
        dfa[j][pat[j]] = j+1
        X = dfa[X][pat[j]]        
    kmp_time += time.clock() - start_time
    return dfa   
    
def combine(str1, str2, overlap):
    return str1 + str2[overlap:]

if __name__ == "__main__":
    main()