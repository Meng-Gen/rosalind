import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    k, n = map(int, lines[0].split())
    array_list = [list(map(int, line.split())) for line in lines[1:]]
    return k, n, array_list

def get_majority_element(array, threshold):
    freq = {}
    for x in array:
        if x not in freq:
            freq[x] = 0
        freq[x] += 1
    for y in freq:
        if freq[y] > threshold:
            return y
    return -1
    
def main():
    k, n, array_list = read_dataset()
    majority_element_list = [get_majority_element(array, n//2) for array in array_list]
    print(' '.join(map(str, majority_element_list)))
    
if __name__ == '__main__':
    sys.exit(main())
