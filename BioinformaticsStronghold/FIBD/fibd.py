import sys
	
def main():
    n, m = map(int, sys.stdin.read().strip().split())
    newborns = [0 for _ in range(m)]
    
    # month #1
    newborns.append(1)
    
    # month #2 to #n
    for i in range(n - 1):
        num_parent = sum([newborns[j] for j in range(-m, -1)])
        newborns.append(num_parent)

    # count alive rabbits
    num_rabbits = sum([newborns[j] for j in range(-m, 0)])
    print(num_rabbits)
        
if __name__ == '__main__':
    sys.exit(main())
