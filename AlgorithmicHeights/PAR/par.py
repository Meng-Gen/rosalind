import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    n = int(lines[0])
    array = list(map(int, lines[1].split()))
    assert(n == len(array))
    return n, array
    
class QuickSort():
    def __init__(self, size, array):
        self.size = size
        self.array = [None] + array
        
    def partition(self):
        self.array[1], self.array[self.size] = self.array[self.size], self.array[1]
        x = self.array[self.size]
        i = 0
        for j in range(1, self.size):
            if self.array[j] <= x:
                i += 1
                self.array[i], self.array[j] = self.array[j], self.array[i]
        self.array[i+1], self.array[self.size] = self.array[self.size], self.array[i+1]
        return self.array[1:]        

    def hoare_partition(self):
        x = self.array[1]
        i, j = 1, self.size
        while True:
            while self.array[j] > x:
                j -= 1
            while self.array[i] < x:
                i += 1
            if i < j:
                self.array[i], self.array[j] = self.array[j], self.array[i]
            else:
                return self.array[1:]
        
def main():
    n, array = read_dataset()
    quick_sort = QuickSort(n, array)
    print(' '.join(map(str, quick_sort.hoare_partition())))
    
if __name__ == '__main__':
    sys.exit(main())
