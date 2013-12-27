import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    k, n = map(int, lines[0].split())
    array_list = [list(map(int, line.split())) for line in lines[1:]]
    return k, n, array_list

def has_two_sum(unsorted_array, t):
    n = len(unsorted_array)
    unsorted_dict = {}
    for i in range(n):
        unsorted_dict[i] = unsorted_array[i]
    sorted_dict = sorted(unsorted_dict.items(), key=lambda x: x[1])
    j, k = 0, n-1
    while j < k:
        sum = sorted_dict[j][1] + sorted_dict[k][1]
        if sum == t:
            index = [sorted_dict[j][0] + 1, sorted_dict[k][0] + 1] # 1-based
            index.sort()
            return index
        elif sum < t:
            j += 1
        else:
            k -= 1
    return [-1]

def main():
    k, n, array_list = read_dataset()
    for array in array_list:
        two_sum = has_two_sum(array, 0)
        print(' '.join(map(str, two_sum)))
    
if __name__ == '__main__':
    sys.exit(main())
