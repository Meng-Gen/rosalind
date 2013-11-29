import sys
	
def main():
	k, m, n = map(int, sys.stdin.read().split())
	total = (k+m+n)*(k+m+n-1)/2
	dominant = k*(k-1)/2 + k*m + k*n + m*(m-1)/2 * 3/4 + m*n * 1/2
	print('%.5f' % float(dominant/total))

if __name__ == '__main__':
    sys.exit(main())
