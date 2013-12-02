import sys
   
def main():
    n = int(sys.stdin.readline().strip())
    # The number of edge of a tree = the number of vertex - 1.
    # Let m be the number of internal vertex, we have 2*(n + m - 1) = n + 3m or m = n - 2.
    print(n - 2)

if __name__ == '__main__':
    sys.exit(main())
