import sys
    
def main():
    n = int(sys.stdin.readline().strip())
    p_array = map(float, sys.stdin.readline().strip().split())
    expected_value_array = [n*p for p in p_array]
    print(' '.join(map(str, expected_value_array)))
    
if __name__ == '__main__':
    sys.exit(main())
