import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    sorted_array = list(map(int, lines[2].split()))
    query_list = list(map(int, lines[3].split()))
    return sorted_array, query_list
    
class BinarySearch():
    def __init__(self, sorted_array):
        self.sorted_array = sorted_array
        self.n = len(sorted_array)

    def find(self, x):
        L, U = 0, self.n - 1
        while L <= U:
            M = (L + U) // 2
            if self.sorted_array[M] < x:
                L = M + 1
            elif self.sorted_array[M] == x:
                return M + 1 # 1-based
            else:
                U = M - 1
        return -1
        
def main():
    sorted_array, query_list = read_dataset()
    search = BinarySearch(sorted_array)
    index_list = [search.find(x) for x in query_list]
    print(' '.join(map(str, index_list)))
        
if __name__ == '__main__':
    sys.exit(main())
