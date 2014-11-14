# Written in Python 3.3.2
# Author:  Rudy Laprade

import heapq, time

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
        print("Time to complete {}:  {} seconds \n".format(i, curr_t-prev_t))
        prev_t = curr_t
    print("Total running time: {} min, {} sec".format(curr_t//60, curr_t%60))

def reconstruct(file):
    dfa_dict = {}
    # def kmp_overlap(str1, str2):
        # alpha = ('A', 'C', 'G', 'T')
        # len2 = len(str2)
        
        # dfa = [{} for _ in range(len2)]
        # for char in alpha: dfa[0][char] = 0
        # dfa[0][str2[0]] = 1
        # X = 0
        
        # j = 0
        # for c in str1:
            # if j == len2:
                # return j
            # if not dfa[j]:
                # for char in alpha:
                    # dfa[j][char] = dfa[X][char]
                # dfa[j][str2[j]] = j+1
                # X = dfa[X][str2[j]]  
            # j = dfa[j][c]
        # dfa_dict[str2] = dfa
        # return j
        
    def kmp_overlap(str1, str2):
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
        return j
    
        
    over_dict = {}
    def overlap(str1, str2):
        ans = kmp_overlap(str1, str2)
        over_dict[(str1, str2)] = ans
        return ans

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
            if i != j:
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
            k = key
            k_len = len(newStr)
            key += 1
            l_str = strand[i]
            r_str = strand[j]
            # Calculate overlaps between merged strand and other strands
            for index in range(k):
                other = strand[index]
                if other != None:
                    heapq.heappush(pq, (-overlap(other,l_str),(index,k)))
                    heapq.heappush(pq, (-overlap(r_str, other),(k,index)))
                    rem = other[overlap(l_str,other):]
                    if rem == r_str[over:over+len(rem)]:
                        heapq.heappush(pq, (-len(other),(k,index)))
            strand[i], strand[j] = None, None # invalidate previous strands
    print("Time to compute other part {}".format(time.clock()-begin))
    return strand[key-1]
    
def kmp(pat):
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
    return dfa   
    
def kmp_overlap_out(str1, str2):
    # dfa = kmp(str2)
    alpha = ('A', 'C', 'G', 'T')
    len2 = len(str2)
    dfa = [{} for _ in range(len2)]
    for char in alpha: dfa[0][char] = 0
    dfa[0][str2[0]] = 1
    X = 0
    
    j = 0
    for c in str1:
        if j == len2:
            return j
        if not dfa[j]:
            for char in alpha:
                dfa[j][char] = dfa[X][char]
            dfa[j][str2[j]] = j+1
            X = dfa[X][str2[j]]  
        j = dfa[j][c]
    return j
    
def combine(str1, str2, overlap):
    return str1 + str2[overlap:]

if __name__ == "__main__":
    main()