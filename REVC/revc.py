import sys

revc = sys.stdin.read().strip()
print(revc.replace("A", "a").replace("T", "A").replace("a", "T").replace("C", "c").replace("G", "C").replace("c", "G")[::-1])
