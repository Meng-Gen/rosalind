import itertools
import sys
    
def minkowski_difference(s, t):
    return [x-y for x, y in itertools.product(s, t)]

def get_frequency(multiset):
    frequency = {}
    for x in multiset:
        entry = None
        for y in frequency:
            if abs(x - y) < 1e-6:
                entry = y
                break
        if not entry:
            frequency[x] = 1
        else:
            frequency[entry] += 1
    return frequency

def get_max_frequency(frequency):
    max_so_far, max_key = 0, None
    for key in frequency:
        if frequency[key] > max_so_far:
            max_so_far, max_key = frequency[key], key
    return max_so_far, max_key
    
def main():
    s = list(map(float, sys.stdin.readline().strip().split()))
    t = list(map(float, sys.stdin.readline().strip().split()))
    difference = minkowski_difference(s, t)
    frequency = get_frequency(difference)
    print('\n'.join(map(str, get_max_frequency(frequency))))
    
if __name__ == '__main__':
    sys.exit(main())
