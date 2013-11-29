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

def first_diff_pos(p, q):
    for i in range(len(p)):
        if p[i] != q[i]:
            return i
            
def reversal_distance(p, q):
    count = 0
    while p != q:
        i = first_diff_pos(p, q)
        j = p.index(q[i])
        p = p[:i] + p[i:j+1][::-1] + p[j+1:]
        count += 1
    return count
    
def main():
    pairs = read_permutation_pairs()
    distances = [reversal_distance(_[0], _[1]) for _ in pairs]
    print(' '.join(map(str, distances)))

if __name__ == '__main__':
    sys.exit(main())
