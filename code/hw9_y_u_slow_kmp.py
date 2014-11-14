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
    over_dict = {}
    def overlap(str1, str2):
        # if (str1, str2) in over_dict:
            # return over_dict[(str1, str2)]
        # ans = 0
        # for i in range(len(str1)):
            # if len(str1)-i > len(str2):
                # if str1[i:i+len(str2)] == str2:
                    # ans = len(str2)
                    # break
            # elif str1[i:] == str2[:len(str1)-i]:
                # ans = len(str1)-i
                # break
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
        # If both strands valid, merge them and add merged strand to end of list
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
            i_len = len(strand[i])
            j_len = len(strand[j])
            # Calculate overlaps between merged strand and other strands
            for index in range(k):
                other = strand[index]
                # if other != None:
                    # heapq.heappush(pq, (-overlap(newStr,other),(k,index)))
                    # heapq.heappush(pq, (-overlap(other,newStr),(index,k)))
                if other != None:
                    heapq.heappush(pq, (-overlap(other,l_str),(index,k)))
                    heapq.heappush(pq, (-overlap(r_str, other),(k,index)))
                    rem = other[overlap(l_str,other):]
                    if rem == r_str[over:over+len(rem)]:
                        heapq.heappush(pq, (-len(other),(k,index)))
                # if other != None and index != i and index != j:
                    # l = len(other)
                    # if overlap(strand[i], other) + overlap(other, strand[j]) >= l:
                        # heapq.heappush(pq, (l,(k, index)))
                    # elif overlap(strand[i], other) + overlap(other, strand[j]) >= k_len:
                        # heapq.heappush(pq, (k_len,(index,k)))
                    # else:
                        # if overlap(other, strand[i]) <= i_len:
                            # heapq.heappush(pq, (-overlap(strand[j],other),(k,index)))
                        # else:
                            # # heapq.heappush(pq, (0,(k,index)))
                            # heapq.heappush(pq, (-overlap(newStr,other),(k,index)))
                        # if overlap(other, strand[j]) <= j_len:
                            # heapq.heappush(pq, (-overlap(other,strand[i]),(index,k)))
                        # else:
                            # # heapq.heappush(pq, (0,(index,k)))
                            # heapq.heappush(pq, (-overlap(other,newStr),(index,k)))
            strand[i], strand[j] = None, None # invalidate previous strands
    print("Time to compute other part {}".format(time.clock()-begin))
    return strand[k]
    
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
    
def kmp_overlap(str1, str2):
    dfa = kmp(str2)
    len2 = len(str2)
    j = 0
    for i in range(len(str1)):
        if j == len2:
            return j
        j = dfa[j][str1[i]]
    return j
    
def combine(str1, str2, overlap):
    return str1 + str2[overlap:]

if __name__ == "__main__":
    main()