import sys
	
def main():
	n, k = map(int, sys.stdin.read().split())
	newborn = [1]
	productive = [0]
	for i in range(n-1):
		next_newborn = k * productive[-1]
		next_productive = productive[-1] + newborn[-1]
		newborn.append(next_newborn)
		productive.append(next_productive)
	print(newborn[-1] + productive[-1])
		
if __name__ == '__main__':
    sys.exit(main())
