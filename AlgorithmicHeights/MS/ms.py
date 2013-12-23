import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    n = int(lines[0])
    unsorted_array = list(map(int, lines[1].split()))
    return n, unsorted_array

def main():
    n, unsorted_array = read_dataset()
    unsorted_array.sort()
    print(' '.join(map(str, unsorted_array)))
    
if __name__ == '__main__':
    sys.exit(main())
