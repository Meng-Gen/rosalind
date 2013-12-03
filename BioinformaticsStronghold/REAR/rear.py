import itertools
import sys

def read_permutation_pairs():
    pairs = []
    is_first_permutation = True
    p, q = None, None
    for line in sys.stdin.readlines():
        line = line.strip()
        if not line:
            continue
        if is_first_permutation:
            p = list(map(int, line.split()))
        else:
            q = list(map(int, line.split()))
            pairs.append([p, q])
        is_first_permutation = not is_first_permutation
    return pairs

def normalize(pair):
    p, q, r = pair[0], pair[1], []
    for i in range(len(q)):
        for j in range(len(p)):
            if p[j] == q[i]:
                r.append(j+1)
                break
    return r

# Reference Book: 
#   AN INTRODUCTION TO BIOINFORMATICS ALGORITHMS, NEIL C. JONES AND PAVEL A. PEVZNER 
def num_breakpoints(permutation):
    num = 0
    prev_c = 0
    for c in permutation:
        if c != prev_c + 1 and c != prev_c - 1:
            num += 1
        prev_c = c
    if prev_c != len(permutation):
        num += 1
    return num

def has_decreasing_strip(permutation):
    return True

def get_min_reversal_permutation(permutation):
    n = len(permutation)
    min_num_breakpoints = n + 1
    min_permutation = None
    for r in itertools.combinations(range(n), 2):
        next = get_reversal_permutation(permutation, r)
        next_num_breakpoints = num_breakpoints(next)
        if min_num_breakpoints > next_num_breakpoints:
            min_num_breakpoints = next_num_breakpoints
            min_permutation = next
    return min_permutation

def flip_increasing_strip(permutation):
    strip = []
    for i in range(1, len(permutation)):
        if permutation[i] + 1 == permutation[i-1]:
            if not strip:
                strip.append(permutation[i-1])
            strip.append(permutation[i])
        else:
            break
    
def get_reversal_permutation(permutation, r):
    return permutation[:r[0]] + permutation[r[0]:r[1]+1][::-1] + permutation[r[1]+1:]

def breakpoint_reversal_distance(permutation):
    n = len(permutation)
    distance = 0
    print(permutation, '<= BEGIN')
    while num_breakpoints(permutation) > 0:
        if has_decreasing_strip(permutation):
            permutation = get_min_reversal_permutation(permutation)
        else:
            #print('''Could not found decreasing strip''')
            #flip_increasing_strip(permutation)
            permutation = get_min_reversal_permutation(permutation)
        print(permutation)
        distance += 1
    return distance
    
def main():
    pairs = read_permutation_pairs()
    normalized_permutations = [normalize(pair) for pair in pairs]
    distances = []
    for p in normalized_permutations:
        distances.append(breakpoint_reversal_distance(p))
    print(' '.join(map(str, distances)))
    
if __name__ == '__main__':
    sys.exit(main())
