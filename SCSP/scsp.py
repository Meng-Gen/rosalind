import itertools
import sys

def shortest_common_supersequence(s, t):
    # Avoid RuntimeError: maximum recursion depth exceeded in comparison
    sys.setrecursionlimit(1500)
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
    if i == 0:
        return t[0:j]
    if j == 0:
        return s[0:i]
    if s[i-1] == t[j-1]:
        return shortest_common_supersequence_printable(c, s, t, i-1, j-1) + s[i-1]
    if c[i][j-1] > c[i-1][j]:
        return shortest_common_supersequence_printable(c, s, t, i, j-1) + t[j-1]
    else:
        return shortest_common_supersequence_printable(c, s, t, i-1, j) + s[i-1]

def main():
    s = sys.stdin.readline().strip()
    t = sys.stdin.readline().strip()
    print(shortest_common_supersequence(s, t))

if __name__ == '__main__':
    sys.exit(main())
