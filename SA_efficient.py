import time
import timeit

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

def efficient_prefix(x, y, simMatrix, gapPenalty):
    n, m = len(x), len(y)
    # mat = [[0 for i in range(2)] for j in range(n + 1)]
    mat = [[0 for i in range(m + 1)] for j in range(2)]
    #initialize the base cases 
    for i in range(m+1):
        mat[0][i] = i * gapPenalty
    #update the matrix in row order 
    for i in range(1, n+1):
        mat[1][0] = mat[0][0] + gapPenalty
        for j in range(1, m+1):
            mat[1][j] = min(mat[0][j-1] + simMatrix[x[i-1] + y[j-1]],
                            mat[0][j] + gapPenalty,
                            mat[1][j-1] + gapPenalty)
        
        for i in range(0, m + 1):
            mat[0][i] = mat[1][i]
    
    return mat[1]

def efficient_suffix(x, y, simMatrix, gapPenalty):
    n, m = len(x), len(y)
    # mat = [[0 for i in range(2)] for j in range(n + 1)]
    mat = [[0 for i in range(m + 1)] for j in range(2)]
    #initialize the base cases 
    for i in range(m+1):
        mat[0][i] = i * gapPenalty
    #update the matrix in row order 
    for i in range(1, n+1):
        mat[1][0] = mat[0][0] + gapPenalty
        for j in range(1, m+1):
            mat[1][j] = min(mat[0][j-1] + simMatrix[x[n-i] + y[m-j]],
                            mat[0][j] + gapPenalty,
                            mat[1][j-1] + gapPenalty)
        
        for i in range(0, m + 1):
            mat[0][i] = mat[1][i]
    
    return mat[1]

def space_efficient_alignment(x, y, simMatrix, gapPenalty):
    # This is the main space_efficient_alignment routine.
    n, m = len(x), len(y)
    if n<2 or m<2:
        # In this case we just use the N-W algorithm.
        # return nw(x, y, simMatrix, gapPenalty, alphEnum)
        return SequenceAlignment(x, y)

    else:
        # Make partitions, call subroutines.
        # F, B = forwards(x[:n//2], y, simMatrix, gapPenalty, alphEnum), backwards(x[n//2:], y, simMatrix, gapPenalty, alphEnum)
        # F, B = forwards(x[:n//2], y, simMatrix, gapPenalty, alphEnum), backwards(x[n//2:], y, simMatrix, gapPenalty, alphEnum)
        F, B = efficient_prefix(x[:n//2], y, simMatrix, gapPenalty), efficient_suffix(x[n//2:], y, simMatrix, gapPenalty)
        # assert(E == F)
        # assert(A == B)
        partition = [F[j] + B[m-j] for j in range(m+1)]
        cut = partition.index(min(partition))
        # Clear all memory now, so that we don't store data during recursive calls.
        F, B, partition = [], [], []
        # Now make recursive calls.
        callLeft = space_efficient_alignment(x[:n//2], y[:cut], simMatrix, gapPenalty)
        callRight = space_efficient_alignment(x[n//2:], y[cut:], simMatrix, gapPenalty)
        # Now return result in format: [1st alignment, 2nd alignment, similarity]
        return [callLeft[r] + callRight[r] for r in range(3)]

def runtime():
    set_up = '''
from __main__ import SequenceAlignment
s1, e1, s2, e2 = read_file('input2.txt')
s1 = expand(s1, e1)
s2 = expand(s2, e2)
gapPenalty = 30
simMatrix = {'AA': 0, 'AC':110, 'CA': 110, 'AG': 48, 'GA': 48, 'CC':0, 'AT': 94, 'TA': 94, 'CG':118, 'GC': 118, 'GG': 0, 'TG':110, 'GT':110, 'TT':0, 'CT': 48, 'TC':48}'''

    test = '''
space_efficient_alignment(s1, s2, simMatrix, gapPenalty) 
'''
    print(timeit.timeit(stmt=test, setup=set_up, globals=globals(), number=1000))

if __name__ == "__main__":
    # s1 = array([G, T, A, C, A, G, T, A], dtype=np.int16)
    # s2 = array([G, G, T, A, C, G, T], dtype=np.int16)
    # aligner = AlignmentFinder(s1, s2)
    # pairs = aligner.find_gobal_alignment()
    # print_sequences(pairs)
    gapPenalty = 30
    simMatrix = {'AA': 0, 'AC':110, 'CA': 110, 'AG': 48, 'GA': 48, 'CC':0, 'AT': 94, 'TA': 94, 'CG':118, 'GC': 118, 'GG': 0, 'TG':110, 'GT':110, 'TT':0, 'CT': 48, 'TC':48}
    s1, e1, s2, e2 = read_file('input2.txt')
    s1 = expand(s1, e1)
    # print(s1)
    s2 = expand(s2, e2)
    # alphEnum = ''
    # s1 = 'ACCT'
    # s2 = 'ACGT'
    z = space_efficient_alignment(s1, s2, simMatrix, gapPenalty)
    print(time.process_time())
    print("Alignment of A: ", z[0])
    print("Alignment of B: ", z[1])
    print("Similarity score: ", z[2], '\n')
    # write_file(z[0], z[1], z[2], runtime, memory_used)
    #/usr/bin/time -l -h python3 SA_efficient.py
    # runtime()