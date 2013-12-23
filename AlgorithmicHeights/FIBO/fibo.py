import sys
    
def fibonacci(n):
    x, y = 0, 1
    for _ in range(n):
        x, y = y, x+y
    return x

def main():
    n = int(sys.stdin.readline().strip())
    print(fibonacci(n))
        
if __name__ == '__main__':
    sys.exit(main())
