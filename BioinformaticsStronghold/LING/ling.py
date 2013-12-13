import math
import sys

def get_substring_count(string):
    n = len(string)
    suffix_list = []
    for i in range(n):
        suffix_list.append(string[i:])
    suffix_list.sort()
    substring_count = n * (n+1) // 2
    for i in range(1, n):
        substring_count -= get_longest_common_prefix(suffix_list[i-1], suffix_list[i])
    return substring_count
    
def get_maximum_substring_count(string):
    n = len(string)
    threshold = int(math.ceil(math.log(n)/math.log(4)))
    return sum([4**k for k in range(1, threshold)]) \
            + sum([n - k + 1 for k in range(threshold, n + 1)])
    
def get_linguistic_complexity(string):
    substring_count = get_substring_count(string)
    #print(substring_count)
    maximum_substring_count = get_maximum_substring_count(string)
    #print(maximum_substring_count)
    return float(substring_count) / float(maximum_substring_count)
        
def get_longest_common_prefix(x, y):
    n = min(len(x), len(y))
    longest_common_pos = n
    for i in range(n):
        if x[i] != y[i]:
            longest_common_pos = i
            break
    return longest_common_pos
    
def main():
    dna_string = sys.stdin.readline().strip()
    print(get_linguistic_complexity(dna_string))
    
if __name__ == '__main__':
    sys.exit(main())
