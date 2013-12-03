import math
import sys

def C(n, m):
    return math.factorial(n) // math.factorial(m) // math.factorial(n - m)

def compute_prob(k, N):
    prob = 0.0
    num_children = 2**k
    for i in range(N, num_children + 1):
        prob += C(num_children, i) * 0.25**i * 0.75**(num_children - i)
    return prob

def main():
    k, N = map(int, sys.stdin.readline().strip().split())
    print(compute_prob(k, N))
    
if __name__ == '__main__':
    sys.exit(main())
