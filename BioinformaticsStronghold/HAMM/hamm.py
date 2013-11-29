import sys

s = sys.stdin.readline().strip()
t = sys.stdin.readline().strip()

print(sum([1 if s[i] != t[i] else 0 for i in range(0, len(s))]))
