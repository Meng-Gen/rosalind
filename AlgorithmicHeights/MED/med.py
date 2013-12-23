import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    n, k = int(lines[0]), int(lines[2])
    unsorted_array = list(map(int, lines[1].split()))
    return n, unsorted_array, k

def main():
    n, unsorted_array, k = read_dataset()
    unsorted_array.sort()
    print(unsorted_array[k-1])
    
if __name__ == '__main__':
    sys.exit(main())
