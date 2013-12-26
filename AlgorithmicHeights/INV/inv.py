import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    n = int(lines[0])
    array = list(map(int, lines[1].split()))
    assert(n == len(array))
    return n, array
    
class CoutingInversionProblem():
    def __init__(self, size, array):
        self.size = size
        self.array = [None] + array
        
    def solve(self):
        return self.count(1, self.size)
    
    def count(self, left, right):
        if left >= right:
            return 0
        middle = (left + right) // 2
        return self.count(left, middle) + self.count(middle + 1, right) + self.merge(left, middle, right)
        
    def merge(self, left, middle, right):
        m = middle - left + 1
        n = right - middle
        L = [None] + self.array[left:middle+1]
        R = [None] + self.array[middle+1:right+1]
        rv = 0
        i, j, k = 1, 1, left
        while i <= m and j <= n:
            if L[i] <= R[j]:
                self.array[k] = L[i]
                i += 1
                rv += (j - 1)
            else:
                self.array[k] = R[j]
                j += 1
            k += 1
        if i > m:
            for y in range(j, n + 1):
                self.array[middle + y] = R[y]
        else:
            for x in range(i, m + 1):
                self.array[right - m + x] = L[x]
            rv += (j - 1)*(m - i + 1)
        return rv
        
def main():
    n, array = read_dataset()
    problem = CoutingInversionProblem(n, array)
    print(problem.solve())
    
if __name__ == '__main__':
    sys.exit(main())
