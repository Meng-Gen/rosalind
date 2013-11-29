import bisect
import sys

def longest_increasing_subsequence(sequence):
    # Robinson-Schensted-Knuth Algorithm
    # Time complexity is O(n log(n))
    if not sequence:
        return []
    aux_len = [1]
    aux_rv = [sequence[0]] 
    for i in range(1, len(sequence)):
        x = sequence[i]
        if x > aux_rv[-1]:
            aux_rv.append(x)
            aux_len.append(len(aux_rv))
        else:
            j = bisect.bisect_right(aux_rv, x)
            aux_rv[j] = x
            aux_len.append(j + 1)        
    rv = []
    curr_len = len(aux_rv)
    for i in range(len(sequence) - 1, -1, -1):
        if curr_len == aux_len[i]:
            rv.append(sequence[i])
            curr_len -= 1
    return rv[::-1]
    
def longest_decreasing_subsequence(sequence):
    return longest_increasing_subsequence(sequence[::-1])[::-1]

def main():
    n = int(sys.stdin.readline().strip())
    p = list(map(int, sys.stdin.readline().strip().split()))
    assert(n == len(p))
    print(' '.join(map(str, longest_increasing_subsequence(p))))
    print(' '.join(map(str, longest_decreasing_subsequence(p))))

if __name__ == '__main__':
    sys.exit(main())