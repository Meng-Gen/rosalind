import sys

def parse():
    pattern = sys.stdin.readline().strip()
    text = sys.stdin.readline().strip()
    d = int(sys.stdin.readline().strip())
    return pattern, text, d

def is_approximating_matching(p, q, d):
    return sum([1 if p[i] != q[i] else 0 for i in range(len(p))]) <= d
    
def main():
    pattern, text, d = parse()
    all_positions = []
    m, n = len(pattern), len(text)
    for i in range(n):
        motif = text[i:i+m]
        if len(motif) == m and is_approximating_matching(motif, pattern, d):
            all_positions.append(i)
    print(' '.join(map(str, all_positions)))

if __name__ == '__main__':
    sys.exit(main())