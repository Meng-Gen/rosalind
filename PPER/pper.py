import math
import sys

def main():
    n, k = map(int, sys.stdin.read().strip().split())
    rv = math.factorial(n) // math.factorial(n-k) % 1000000
    print(rv)
    
if __name__ == '__main__':
    sys.exit(main())
