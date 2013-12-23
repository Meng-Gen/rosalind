import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    array_list = [list(map(int, line.split())) for line in lines]
    return array_list[1], array_list[3]

def merge(A, B):
    m, n = len(A), len(B)
    i, j = 0, 0
    C = []
    while i < m and j < n:
        if A[i] < B[j]:
            C.append(A[i])
            i += 1
        else:
            C.append(B[j])
            j += 1
    if i == m:
        C += B[j:]
    if j == n:
        C += A[i:]
    return C
    
def main():
    A, B = read_dataset()
    C = merge(A, B)
    print(' '.join(map(str, C)))
    
if __name__ == '__main__':
    sys.exit(main())
