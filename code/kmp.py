
alpha = ('A', 'C', 'G', 'T')


def kmp(pat):
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
    for c in str1:
        if j == len2:
            return j
        j = dfa[j][c]
    return j