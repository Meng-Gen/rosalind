import itertools
import sys

def shortest_common_supersequence(s, t):
    c = longest_common_subsequence_table(s, t)
    return shortest_common_supersequence_printable(c, s, t, len(s), len(t))

def longest_common_subsequence_table(s, t):
    m, n = len(s), len(t)
    c = [[0 for _ in range(n+1)] for _ in range(m+1)] 
    for i in range(m):
        for j in range(n):
            if s[i] == t[j]:
                c[i+1][j+1] = c[i][j] + 1
            else:
                c[i+1][j+1] = max(c[i+1][j], c[i][j+1])
    return c

def shortest_common_supersequence_printable(c, s, t, i, j):
    reversed_scs = ''
    while True:
        if i == 0:
            reversed_scs += t[0:j][::-1]
            break
        if j == 0:
            reversed_scs += s[0:i][::-1]
            break
        if s[i-1] == t[j-1]:
            reversed_scs += s[i-1]
            i, j = i-1, j-1
        elif c[i][j-1] > c[i-1][j]:
            reversed_scs += t[j-1]
            i, j = i, j-1
        else:
            reversed_scs += s[i-1]
            i, j = i-1, j
    return reversed_scs[::-1]

def main():
    s = sys.stdin.readline().strip()
    t = sys.stdin.readline().strip()
    print(shortest_common_supersequence(s, t))

if __name__ == '__main__':
    sys.exit(main())
