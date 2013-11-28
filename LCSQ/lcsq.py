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
    c = longest_common_subsequence_table(s, t)
    return longest_common_subsequence_backtrack(c, s, t, len(s), len(t))

def longest_common_subsequence_table(s, t):
    m, n = len(s), len(t)
    c = [[0 for x in range(n+1)] for y in range(m+1)] 
    for i in range(m):
        for j in range(n):
            if s[i] == t[j]:
                c[i+1][j+1] = c[i][j] + 1
            else:
                c[i+1][j+1] = max(c[i+1][j], c[i][j+1])
    return c

# May introduce RuntimeError: maximum recursion depth exceeded in comparison
def longest_common_subsequence_backtrack(c, s, t, i, j):
    if i == 0 or j == 0:
        return ''
    if s[i-1] == t[j-1]:
        return longest_common_subsequence_backtrack(c, s, t, i-1, j-1) + s[i-1]
    if c[i][j-1] > c[i-1][j]:
        return longest_common_subsequence_backtrack(c, s, t, i, j-1)
    else:
        return longest_common_subsequence_backtrack(c, s, t, i-1, j)
            
def main():
    # Avoid RuntimeError: maximum recursion depth exceeded in comparison
    sys.setrecursionlimit(1500)
    fasta_records = get_fasta_records(sys.stdin.readlines())
    s, t = fasta_records[0][1], fasta_records[1][1]
    print(longest_common_subsequence(s, t))

if __name__ == '__main__':
    sys.exit(main())
