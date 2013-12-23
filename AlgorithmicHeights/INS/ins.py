import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    return list(map(int, lines[1].split()))
    
class InsertionSorter():
    def __init__(self, unsorted_array):
        self.unsorted_array = unsorted_array
        self.sorted_array = None
        self.num_swaps = None
    
    def get_num_swaps(self):
        self.sort()
        return self.num_swaps
        
    def sort(self):
        self.num_swaps = 0
        self.sorted_array = list(self.unsorted_array)
        n = len(self.sorted_array)
        for i in range(1, n):
            k = i
            while k > 0 and self.sorted_array[k] < self.sorted_array[k-1]:
                self.sorted_array[k], self.sorted_array[k-1] = self.sorted_array[k-1], self.sorted_array[k]
                k -= 1
                self.num_swaps += 1
                
def main():
    array = read_dataset()
    sorter = InsertionSorter(array)
    print(sorter.get_num_swaps())
    
if __name__ == '__main__':
    sys.exit(main())
