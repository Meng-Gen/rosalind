import sys

def b(n):
    count, num_edge = 1, 3
    for i in range(3, n):
        count = (count * num_edge) % 1000000
        num_edge += 2
    return count

def main():
    n = int(sys.stdin.readline().strip())
    print(b(n))
    
if __name__ == '__main__':
    sys.exit(main())
