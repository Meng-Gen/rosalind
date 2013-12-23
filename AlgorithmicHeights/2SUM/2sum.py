import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    k, n = map(int, lines[0].split())
    array_list = [list(map(int, line.split())) for line in lines[1:]]
    return k, n, array_list

def has_two_sum(unsorted_array, t):
    n = len(unsorted_array)
    unsorted_array_map = {}
    for i in range(n):
        x = unsorted_array[i] 
        if x not in unsorted_array_map:
            unsorted_array_map[x] = set()
        unsorted_array_map[x].add(i)
    A = list(unsorted_array)
    A.sort()
    j, k = 0, len(A)-1
    while j < k:
        sum = A[j] + A[k]
        if sum == t:
            p, q = None, None
            if A[j] * 2 == sum:
                p, q = list(unsorted_array_map[A[j]])[:2]
            else:
                p = list(unsorted_array_map[A[j]])[0]
                q = list(unsorted_array_map[A[k]])[0]
            if p < q:
                return [p+1, q+1] # 1-based
            else:
                return [q+1, p+1] # 1-based
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
