import sys

a, b = map(int, sys.stdin.readline().split())
print(sum(_ if _ % 2 == 1 else 0 for _ in range(a, b+1)))
