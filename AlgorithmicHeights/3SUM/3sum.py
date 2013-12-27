import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    k, n = map(int, lines[0].split())
    array_list = [list(map(int, line.split())) for line in lines[1:]]
    return k, n, array_list

def has_three_sum(unsorted_array):
    n = len(unsorted_array)
    unsorted_dict = {}
    for i in range(n):
        unsorted_dict[i] = unsorted_array[i]
    sorted_dict = sorted(unsorted_dict.items(), key=lambda x: x[1])
    
    for i in range(n):
        j, k = i+1, n-1
        t = -sorted_dict[i][1]
        while j < k:
            sum = sorted_dict[j][1] + sorted_dict[k][1]
            if sum == t:
                index = [sorted_dict[i][0] + 1, sorted_dict[j][0] + 1, sorted_dict[k][0] + 1] # 1-based
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
        three_sum = has_three_sum(array)
        print(' '.join(map(str, three_sum)))
    
if __name__ == '__main__':
    sys.exit(main())
