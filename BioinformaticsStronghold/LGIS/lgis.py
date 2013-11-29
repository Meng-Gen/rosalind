import sys

def LIS(sequence):
    # Time complexity is O(n^2)	
    lis_len = [1] * len(sequence)
    for i in range(len(sequence)):
        for j in range(i+1, len(sequence)):
            if sequence[i] < sequence[j]:
                lis_len[j] = max(lis_len[j], lis_len[i] + 1)
    rv = []
    curr_len = max(lis_len)
    for i in range(len(sequence) - 1, -1, -1):
        if curr_len == lis_len[i]:
            rv.append(sequence[i])
            curr_len -= 1
    return rv[::-1]

def LDS(sequence):
    return LIS(sequence[::-1])[::-1]
    
def main():
    n = int(sys.stdin.readline().strip())
    p = list(map(int, sys.stdin.readline().strip().split()))
    assert(len(p) == n)
    print(' '.join(map(str, LIS(p))))
    print(' '.join(map(str, LDS(p))))
		
if __name__ == '__main__':
    sys.exit(main())
