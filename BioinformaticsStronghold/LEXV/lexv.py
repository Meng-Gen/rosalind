import itertools
import sys

def normalize(s):
    head = ''.join(itertools.takewhile(lambda x: x != '!', s))
    tail = ''.join(itertools.takewhile(lambda x: x == '!', reversed(s)))
    if len(head) + len(tail) == len(s):
        return ''.join(head)

def main():
    alphabets = ['!'] + sys.stdin.readline().strip().split()
    n = int(sys.stdin.readline().strip())
    for s in itertools.product(alphabets, repeat=n):
        good_permutation = normalize(s)
        if good_permutation:
            print(good_permutation)
    
if __name__ == '__main__':
    sys.exit(main())
