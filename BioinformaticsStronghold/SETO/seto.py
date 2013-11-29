import sys
   
def main():
    n = int(sys.stdin.readline().strip())
    U = set(range(1, n+1))
    A = eval(sys.stdin.readline().strip())
    B = eval(sys.stdin.readline().strip())
    print(A | B)
    print(A & B)
    print(A - B)
    print(B - A)
    print(U - A)
    print(U - B)

if __name__ == '__main__':
    sys.exit(main())
