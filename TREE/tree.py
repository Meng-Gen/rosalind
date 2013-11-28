import sys
    
def main():
    lines = sys.stdin.readlines()
    n = int(lines[0])
    print(n - len(lines))

if __name__ == '__main__':
    sys.exit(main())
