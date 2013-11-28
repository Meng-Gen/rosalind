import sys

def get_fasta_records(lines):
    dna_strings = []
    is_first_dna = True
    id, dna = '', ''
    for line in lines:
        line = line.strip();
        if line[0] == '>':
            if is_first_dna:
                is_first_dna = False
            else:
                dna_strings.append([id, dna])
            id, dna = line[1:], ''
        else:
            dna += line
    dna_strings.append([id, dna])
    return dna_strings

def longest_common_subsequence(s, t):
    # Avoid RuntimeError: maximum recursion depth exceeded in comparison
    sys.setrecursionlimit(1500)
    c, b = longest_common_subsequence_table(s, t)
    return longest_common_subsequence_printable(b, s, len(s), len(t))

def longest_common_subsequence_table(s, t):
    m, n = len(s), len(t)
    b = [[0 for _ in range(n)] for _ in range(m)] 
    c = [[0 for _ in range(n+1)] for _ in range(m+1)] 
    for i in range(m):
        for j in range(n):
            if s[i] == t[j]:
                c[i+1][j+1] = c[i][j] + 1
                b[i][j] = 'Diagonal'
            elif c[i][j+1] > c[i+1][j]:
                c[i+1][j+1] = c[i][j+1]
                b[i][j] = 'Up'
            else:
                c[i+1][j+1] = c[i+1][j]
                b[i][j] = 'Left'
    return c, b

def longest_common_subsequence_printable(b, s, i, j):
    if i == 0 or j == 0:
        return ''
    direction = b[i-1][j-1]
    if direction == 'Diagonal':
        return longest_common_subsequence_printable(b, s, i-1, j-1) + s[i-1]
    elif direction == 'Up':
        return longest_common_subsequence_printable(b, s, i-1, j)
    else:
        return longest_common_subsequence_printable(b, s, i, j-1)
            
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    s, t = fasta_records[0][1], fasta_records[1][1]
    print(longest_common_subsequence(s, t))

if __name__ == '__main__':
    sys.exit(main())
