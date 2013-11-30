import math
import sys

def get_common_log(content_count, x):
    content_probability = [(1-x)/2.0, x/2.0, x/2.0, (1-x)/2.0]
    return sum([_[0] * math.log10(_[1]) for _ in zip(content_count, content_probability)])

def main():
    dna = sys.stdin.readline().strip()
    probabilities = map(float, sys.stdin.readline().strip().split())
    content_count = [dna.count(i) for i in 'ACGT']    
    all_common_logs = []
    for x in probabilities:
        all_common_logs.append(get_common_log(content_count, x))
    print(' '.join(map(str, all_common_logs)))
        
if __name__ == '__main__':
    sys.exit(main())
