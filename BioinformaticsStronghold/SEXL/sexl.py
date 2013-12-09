import sys
    
def main():
    dataset = map(float, sys.stdin.readline().strip().split())
    prob_array = [2.0 * x * (1.0 - x) for x in dataset]
    print(' '.join(map(str, prob_array)))
    
if __name__ == '__main__':
    sys.exit(main())
