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
        if (str1, str2) in over_dict:
            return over_dict[(str1, str2)]
        ans = 0
        for i in range(len(str1)):
            if len(str1)-i > len(str2):
                if str1[i:i+len(str2)] == str2:
                    ans = len(str2)
                    break
            elif str1[i:] == str2[:len(str1)-i]:
                ans = len(str1)-i
                break
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
    return strand[k]
    

    
def combine(str1, str2, overlap):
    return str1 + str2[overlap:]

if __name__ == "__main__":
    main()