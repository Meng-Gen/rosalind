import sys

def get_combinations(n):
    curr = [1]
    for i in range(1, n+1):
        next = [1]
        for j in range(1, i):
            next.append((curr[j-1] + curr[j]) % 1000000)
        next.append(1)
        curr = next
    return curr
    
def main():
    n, m = map(int, sys.stdin.readline().strip().split())
    combinations = get_combinations(n)
    print(sum(combinations[m:]) % 1000000)
            
if __name__ == '__main__':
    sys.exit(main())
