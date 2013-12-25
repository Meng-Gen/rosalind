import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    n = int(lines[0])
    unsorted_array = list(map(int, lines[1].split()))
    return n, unsorted_array
    
class Heap():
    def __init__(self, size, array):
        self.size = size
        self.array = [None] + array
        
    def heapsort(self):
        self.build()
        for i in range(self.size, 0, -1):
            self.array[1], self.array[i] = self.array[i], self.array[1]
            self.size -= 1
            self.heapify(1)
        return self.array[1:]
    
    def build(self):
        for i in range(self.size // 2, 0, -1):
            self.heapify(i)
        return self.array[1:]
        
    def heapify(self, i):
        l = self.left(i)
        r = self.right(i)
        largest = None
        if l <= self.size and self.array[l] > self.array[i]:
            largest = l
        else:
            largest = i
        if r <= self.size and self.array[r] > self.array[largest]:
            largest = r
        if largest != i:
            self.array[i], self.array[largest] = self.array[largest], self.array[i]
            self.heapify(largest)
        
    def left(self, i):
        return 2*i
    
    def right(self, i):
        return 2*i + 1    
    
def main():
    n, array = read_dataset()
    heap = Heap(n, array)
    print(' '.join(map(str, heap.heapsort())))
    
if __name__ == '__main__':
    sys.exit(main())
