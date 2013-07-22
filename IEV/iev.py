import sys
	
def main():
	expected_values = [ 1, 1, 1, 0.75, 0.5, 0 ]
	pairings = list(map(int, sys.stdin.read().split()))
	print(sum([pairings[i] * expected_values[i] for i in range(len(expected_values))]) * 2)
	
if __name__ == '__main__':
    sys.exit(main())
