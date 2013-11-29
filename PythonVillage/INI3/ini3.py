import sys

s = sys.stdin.readline().strip()
a, b, c, d = map(int, sys.stdin.readline().split())
print(s[a:b+1], s[c:d+1])
