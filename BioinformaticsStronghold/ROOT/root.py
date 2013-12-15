import sys

def B(n):
    count, num_edge = 3, 5
    for i in range(3, n):
        count = (count * num_edge) % 1000000
        num_edge += 2
    return count

def main():
    n = int(sys.stdin.readline().strip())
    print(B(n))
    
if __name__ == '__main__':
    sys.exit(main())
