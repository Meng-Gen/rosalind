import sys

def get_prob(dna, x, N):
    content_count = [dna.count(i) for i in 'ACGT']   
    content_probability = [(1-x)/2.0, x/2.0, x/2.0, (1-x)/2.0]
    prob = 1.0
    for content in zip(content_count, content_probability):
        prob *= content[1]**content[0]
    return 1.0 - (1.0 - prob)**N

def main():
    N, x = sys.stdin.readline().strip().split()
    N, x = int(N), float(x)
    dna = sys.stdin.readline().strip()
    prob = get_prob(dna, x, N)
    print(prob)
    
if __name__ == '__main__':
    sys.exit(main())
