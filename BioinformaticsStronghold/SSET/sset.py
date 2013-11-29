import sys
    
def main():
    n = int(sys.stdin.read().strip())
    print(2**n % 1000000)

if __name__ == '__main__':
    sys.exit(main())
