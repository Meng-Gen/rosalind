import itertools
import sys

n = sys.stdin.read().strip()
all_perm = list(itertools.permutations(range(1, int(n)+1)))
print(len(all_perm))
for perm in all_perm:
	print(' '.join([str(i) for i in perm]))
