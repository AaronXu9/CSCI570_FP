def read_file(filename):
    with open(filename, 'r') as f:
        s1 = f.readline().strip()
        s2 = ''
        next_lines = f.readlines()
        
        s1_expand = []
        idx = 0
        for line in next_lines:
            if line.strip().isdigit():
                s1_expand.append(int(line))
                idx += 1
            else: 
                s2 = line.strip()
                idx += 1
                break
        
        s2_expand = []
        for i in range(idx, len(next_lines)):
            s2_expand.append(int(next_lines[i]))
    
    return s1, s1_expand, s2, s2_expand

def write_file(a1, a2, cost, runtime, memory_used):

    return 
def expand(s:str, s_expand:list):
    for idx in s_expand:
        s = s[:idx+1] + s + s[idx+1:]
    return s

def SequenceAlignment(s1, s2):
    m = len(s1)
    n = len(s2)

    SA = [[0 for i in range(n + 1)] for j in range(m + 1)]
    delta = 30
    alpha = {'AA': 0, 'AC':110, 'CA': 110, 'AG': 48, 'GA': 48, 'CC':0, 'AT': 94, 'TA': 94, 'CG':118, 'GC': 118, 'GG': 0, 'TG':110, 'GT':110, 'TT':0, 'CT': 48, 'TC':48}
    s1_aligned = ''
    s2_aligned = ''

    #initialization A[i, 0]= iδ, 
    for i in range(m + 1):
        SA[i][0] = i * 30
    
    #initialization A[0, j]= jδ, 
    for j in range(n + 1):
        SA[0][j] = j * 30
    
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            SA[i][j] = min(alpha[s1[i-1]+s2[j-1]] + SA[i-1][j-1], delta + SA[i-1][j], delta + SA[i][j-1])

    i = m
    j = n
    while i and j:
        if SA[i][j] == alpha[s1[i-1]+s2[j-1]] + SA[i-1][j-1]:
            s1_aligned = s1[i-1] + s1_aligned
            s2_aligned = s2[j-1] + s2_aligned
            i -= 1
            j -= 1
        elif SA[i][j] == delta + SA[i-1][j]:
            s1_aligned = s1[i-1] + s1_aligned
            s2_aligned = '_' + s2_aligned
            i -= 1

        elif SA[i][j] == delta + SA[i][j-1]:
            s1_aligned = '_' + s1_aligned
            s2_aligned = s2[j-1] + s2_aligned
            j -= 1
    
    while i:
        s1_aligned = s1[i-1] + s1_aligned
        s2_aligned = '_' + s2_aligned
        i -= 1
    while j:
        s1_aligned = '_' + s1_aligned
        s2_aligned = s2[j-1] + s2_aligned
        j -= 1
    
    # alignment,_ = get_aligned(s1, s2, SA, delta, alpha)
    # return SA[m][n], alignment
    return [s1_aligned, s2_aligned, SA[m][n]]

if __name__ == "__main__":
    # s1 = array([G, T, A, C, A, G, T, A], dtype=np.int16)
    # s2 = array([G, G, T, A, C, G, T], dtype=np.int16)
    # aligner = AlignmentFinder(s1, s2)
    # pairs = aligner.find_gobal_alignment()
    # print_sequences(pairs)
    # gapPenalty = 30
    # simMatrix = {'AA': 0, 'AC':110, 'CA': 110, 'AG': 48, 'GA': 48, 'CC':0, 'AT': 94, 'TA': 94, 'CG':118, 'GC': 118, 'GG': 0, 'TG':110, 'GT':110, 'TT':0, 'CT': 48, 'TC':48}
    s1, e1, s2, e2 = read_file('input2.txt')
    s1 = expand(s1, e1)
    # print(s1)
    s2 = expand(s2, e2)
    alphEnum = ''
    # s1 = 'ACCT'
    # s2 = 'ACGT'
    z = SequenceAlignment(s1, s2)
    print("Alignment of A: ", z[0])
    print("Alignment of B: ", z[1])
    print("Similarity score: ", z[2], '\n')
    write_file(z[0], z[1], z[2], runtime, memory_used)