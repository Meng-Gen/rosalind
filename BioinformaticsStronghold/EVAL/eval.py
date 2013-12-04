import sys

def get_prob(dna, x):
    content_count = [dna.count(i) for i in 'ACGT']   
    content_probability = [(1-x)/2.0, x/2.0, x/2.0, (1-x)/2.0]
    prob = 1.0
    for content in zip(content_count, content_probability):
        prob *= content[1]**content[0]
    return prob

def get_expected_number(n, dna, x):
    prob = get_prob(dna, x)
    return prob * (n - len(dna) + 1)

def main():
    n = int(sys.stdin.readline().strip())
    dna = sys.stdin.readline().strip()
    x_list = list(map(float, sys.stdin.readline().strip().split()))
    expected_number_list = [get_expected_number(n, dna, x) for x in x_list]
    print(' '.join(map(str, expected_number_list)))
    
if __name__ == '__main__':
    sys.exit(main())
