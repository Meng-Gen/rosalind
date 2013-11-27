import itertools
import sys

def main():
    alphabets = ''.join(sys.stdin.readline().strip().split())
    n = int(sys.stdin.readline().strip())
    for s in itertools.product(alphabets, repeat=n):
        print(''.join(s))
    
if __name__ == '__main__':
    sys.exit(main())
