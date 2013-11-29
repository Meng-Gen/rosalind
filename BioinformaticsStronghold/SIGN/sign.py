import itertools
import math
import sys

def main():
    n = int(sys.stdin.readline().strip())
    total_count = 2**n * math.factorial(n)
    print(total_count)
    for i in itertools.permutations([str(_) for _ in range(1, n + 1)], n):
        for j in itertools.product('-+', repeat=n):
            print(' '.join(str(int(''.join(k))) for k in zip(j, i)))
    
if __name__ == '__main__':
    sys.exit(main())
