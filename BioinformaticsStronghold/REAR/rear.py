import itertools
import sys

reversal_distance_table = {}

def permutation_hash(p):
    return '_'.join(map(str, p))

def init_reversal_distance_table():
    #root = list(range(1, 11)))
    
    count = 0
    for it in itertools.permutations(range(1, 11)):
        reversal_distance_table[permutation_hash(it)] = None
        count += 1
        if count % 100000 == 0:
            print(permutation_hash(it))
        #break
    #reversal_distance_table[permutation_hash(root)] = 0
    pass

def read_permutation_pairs():
    pairs = []
    is_first_permutation = True
    p, q = None, None
    for line in sys.stdin.readlines():
        line = line.strip()
        if not line:
            continue
        if is_first_permutation:
            p = list(line.split())
        else:
            q = list(line.split())
            pairs.append([p, q])
        is_first_permutation = not is_first_permutation
    return pairs
            
def reversal_distance(p, q):
    pass
    
def main():
    init_reversal_distance_table()
    pairs = read_permutation_pairs()
    distances = [reversal_distance(_[0], _[1]) for _ in pairs]
    print(' '.join(map(str, distances)))
    #print(reversal_distance_table)
    
if __name__ == '__main__':
    sys.exit(main())
