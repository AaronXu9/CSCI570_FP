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

def expand(s:str, s_expand:list):
    for idx in s_expand:
        s = s[:idx+1] + s + s[idx+1:]
    return s

def get_aligned_pair(s1, s2, i, j):
    # aligned = s[idx-1] if idx > 0 else '_'
    n1 = s1[i-1] if i>0 else '_'
    n2 = s2[j-1] if j>0 else '_'
    return (n1, n2)
    # return aligned

def get_aligned(s1, s2, SA, delta, alpha):
    alignment= []
    i = len(s1)
    j = len(s2)

    while i >0 and j>0:
        if SA[i][j] == alpha[s1[i-1]+s2[j-1]] + SA[i-1][j-1]:
            alignment.append(get_aligned_pair(s1, s2, i, j))
            i -= 1
            j -= 1
        elif SA[i][j] == delta + SA[i-1][j]:
            alignment.append(get_aligned_pair(s1, s2, i, 0))
            i -= 1
        else:
            alignment.append(get_aligned_pair(s1, s2, 0, j))
            j -= 1
    
    while i > 0:
        alignment.append(get_aligned_pair(s1, s2, i, 0))
        i -= 1
    while j > 0:
        alignment.append(get_aligned_pair(s1, s2, 0, j))
        j -= 1

    alignment.reverse()
    return alignment

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
    while i >0 and j >0:
        if SA[i][j] == alpha[s1[i-1]+s2[j-1]] + SA[i-1][j-1]:
            s1_aligned += s1[i-1]
            s2_aligned += s2[j-1]
            i -= 1
            j -= 1
        elif SA[i][j] == delta + SA[i-1][j]:
            s2_aligned += '_'
            s1_aligned += s1[i-1]
            i -= 1
            j -= 1
        elif SA[i][j] == delta + SA[i][j-1]:
            s1_aligned += '_'
            s2_aligned += s2[j-1]
            i -= 1
            j -= 1
    
    while i > 0:
        s1_aligned += s1[i-1]
        s2_aligned += '_'
        i -= 1
    while j > 0:
        s1_aligned += '_'
        s2_aligned += s2[j-1]
        j -= 1
    
    alignment = get_aligned(s1, s2, SA, delta, alpha)
    return SA[m][n], alignment
    # return SA[m][n], s1_aligned[::-1], s2_aligned[::-1]

def SequenceAlignmentBackward(s1, s2):
    m = len(s1)
    n = len(s2)

    SA = [[0 for i in range(n + 1)] for j in range(m + 1)]
    delta = 30
    alpha = {'AA': 0, 'AC':110, 'CA': 110, 'AG': 48, 'GA': 48, 'CC':0, 'AT': 94, 'TA': 94, 'CG':118, 'GC': 118, 'GG': 0, 'TG':110, 'GT':110, 'TT':0, 'CT': 48, 'TC':48}
    aligned = ''

    #initialization A[i, 0]= iδ, 
    for i in range(m + 1):
        SA[m-i][n] = i * 30
    
    #initialization A[0, j]= jδ, 
    for j in range(n + 1):
        SA[m][n-j] = j * 30
    
    for j in range(n-1, -1, -1):
        for i in range(m-1, -1, -1):
            SA[i][j] = min(alpha[s1[i]+s2[j]] + SA[i+1][j+1], delta + SA[i+1][j], delta + SA[i][j+1])

    return SA[0][0]

def Space_Efficient_Alignment(s1, s2):
    m = len(s1)
    n = len(s2)

    SA = [[0 for i in range(2)] for j in range(m + 1)]
    delta = 30
    alpha = {'AA': 0, 'AC':110, 'CA': 110, 'AG': 48, 'GA': 48, 'CC':0, 'AT': 94, 'TA': 94, 'CG':118, 'GC': 118, 'GG': 0, 'TG':110, 'GT':110, 'TT':0, 'CT': 48, 'TC':48}
    aligned = ''

    #initialization A[i, 0]= iδ, 
    for i in range(m + 1):
        SA[i][0] = i * 30
    
    #initialization A[0, j]= jδ, 
    # for j in range(n + 1):
    #     SA[0][j] = j * 30
    
    for j in range(1, n + 1):
        SA[0][1] = j * delta
        for i in range(1, m + 1):
            SA[i][1] = min(alpha[s1[i-1]+s2[j-1]] + SA[i-1][0], delta + SA[i-1][1], delta + SA[i][0])
        
        #Move column 1 of B to column 0 to make room for next iteration: Update B[i, 0]= B[i, 1] for each i
        for i in range(0, m + 1):
            SA[i][0] = SA[i][1]
    
    return [row[0] for row in SA]
    return SA[m][0]

def Space_Efficient_Alignment_Backward(s1, s2):
    m = len(s1)
    n = len(s2)

    SA = [[0 for i in range(2)] for j in range(m + 1)]
    delta = 30
    alpha = {'AA': 0, 'AC':110, 'CA': 110, 'AG': 48, 'GA': 48, 'CC':0, 'AT': 94, 'TA': 94, 'CG':118, 'GC': 118, 'GG': 0, 'TG':110, 'GT':110, 'TT':0, 'CT': 48, 'TC':48}
    aligned = ''

    #initialization A[i, 0]= iδ, 
    for i in range(m + 1):
        SA[m-i][1] = i * 30
    
    for j in range(n-1, -1, -1):
        SA[m][0] = j * delta
        for i in range(m-1, -1, -1):
            SA[i][0] = min(alpha[s1[i]+s2[j]] + SA[i+1][1], delta + SA[i+1][0], delta + SA[i][1])
        
        for i in range(0, m + 1):
            SA[i][1] = SA[i][0]

    return [row[0] for row in SA]
    # return SA[0][0]

def Align(s1, s2):
    m = len(s1)
    n = len(s2)
    global opt_path

    if m <= 2 or n <= 2:
        return SequenceAlignment(s1, s2)
    
    x_prefix = Space_Efficient_Alignment(s1, s2[:n//2])
    x_suffix = Space_Efficient_Alignment_Backward(s1, s2[n//2+1:n+1])

    best = float('inf')
    for q in range(m):
        cost = x_prefix[q] + x_suffix[q+1]
        if cost < best:
                bestq = q
                best = cost
    
    opt_path.append((n//2, bestq))
    Align(s1[:bestq+1], s2[:n//2])
    Align(s1[bestq+1:n+1], s2[n//2+1:n])

    return opt_path

def print_seq(pairs):
    top_seq = ''
    bottom_seq = ''
    for (b, t) in pairs:
        top_seq += t
        bottom_seq += b
    
    print(len(top_seq), top_seq)
    print(len(bottom_seq), bottom_seq)
    return top_seq, bottom_seq

def compute_aglignment_cost(s1, s2, delta, alpha):
    cost = 0
    for i in range(len(s1)):
        # print(s1[i], s2[i])
        if s1[i] == '_' or s2[i] == '_':
            cost += delta
        else:
            cost += alpha[s1[i] + s2[i]]
    return cost

opt_path = []

if __name__ == '__main__':
    delta = 30
    alpha = {'AA': 0, 'AC':110, 'CA': 110, 'AG': 48, 'GA': 48, 'CC':0, 'AT': 94, 'TA': 94, 'CG':118, 'GC': 118, 'GG': 0, 'TG':110, 'GT':110, 'TT':0, 'CT': 48, 'TC':48}
    s1_aligned = ''
    s1, e1, s2, e2 = read_file('BaseTestcases_CS570FinalProject/input2.txt')
    s1 = expand(s1, e1)
    s2 = expand(s2, e2)
    # assert(expand(s1, e1) == 'ACACTGACTACTGACTGGTGACTACTGACTGG')
    # assert(expand(s2, e2) == 'TATTATACGCTATTATACGCGACGCGGACGCG')
    test_seq1 = 'ACCT'
    test_seq2 = 'ACGT'
    # print(SequenceAlignment('ACCT', 'ACGT'))
    cost, alignment = SequenceAlignment(test_seq1,test_seq2)
    # cost, alignment = SequenceAlignment(s1,s2)
    a1, a2 = print_seq(alignment)
    print('aligned cost', compute_aglignment_cost(a1, a2, delta, alpha))
    # cost, s1_aligned, s2_aligned = SequenceAlignment(s1, s2)
    # assert(s1_aligned[:50] == 'ACACACTGAC_TACTGACTGGTG_ACTACTG_ACT_G_GACTGAC_TACT')
    # assert(s1_aligned[-50:] == 'TAT_TA_TTATACG_CTA_TTATACG_CGACGCGGACGC_G_T_ATACG_')

    # assert(s2_aligned[:50] == ) 
    # cost, s1_aligned, s2_aligned = SequenceAlignment(test_seq1, test_seq2)
    # print(SequenceAlignment(s1, s2))
    print(cost)
    # print(s1_aligned)
    # print(s2_aligned)

    # print(Space_Efficient_Alignment(s1, s2))
    # print(SequenceAlignmentBackward(s1, s2))
    # print(Space_Efficient_Alignment_Backward(s1, s2))
    s1_aligned = '' 
    s2_aligned = ''
    
    print('--------------------------------')
    Align(test_seq1, test_seq2)

    print(opt_path)