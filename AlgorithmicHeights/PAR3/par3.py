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
        k = i
        for j in range(1, self.size):
            if self.array[j] < x:
                k += 1
                self.array[j], self.array[k] = self.array[k], self.array[j]
                i += 1
                self.array[i], self.array[k] = self.array[k], self.array[i]
            elif self.array[j] <= x:
                k += 1
                self.array[j], self.array[k] = self.array[k], self.array[j]
        self.array[k+1], self.array[self.size] = self.array[self.size], self.array[k+1]
        return self.array[1:]        

def main():
    n, array = read_dataset()
    quick_sort = QuickSort(n, array)
    print(' '.join(map(str, quick_sort.partition())))
    
if __name__ == '__main__':
    sys.exit(main())
