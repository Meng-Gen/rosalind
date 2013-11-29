import sys

is_even = False
for line in sys.stdin.readlines():
    if not is_even:
        is_even = True
    else:
        line = line.strip()
        print(line)
        is_even = False
