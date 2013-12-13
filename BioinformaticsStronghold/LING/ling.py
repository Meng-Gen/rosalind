import math
import sys

def get_substring_count(string):
    n = len(string)
    substring_count = n * (n+1) // 2
    for alphabet_bucket in 'ACGT':
        count = get_substring_count_helper(string, alphabet_bucket)
        substring_count -= count
    return substring_count

# Avoid out of memory by bucket sort    
def get_substring_count_helper(string, alphabet_bucket):
    suffix_list = []
    for i in range(len(string)):
        suffix = string[i:]
        if suffix[0] == alphabet_bucket:
            suffix_list.append(suffix)
    suffix_list.sort()
    count = 0
    for i in range(1, len(suffix_list)):
        count += get_longest_common_prefix(suffix_list[i-1], suffix_list[i])
    return count
    
def get_maximum_substring_count(string):
    n = len(string)
    threshold = int(math.ceil(math.log(n)/math.log(4)))
    return sum([4**k for k in range(1, threshold)]) \
            + sum([n - k + 1 for k in range(threshold, n + 1)])
    
def get_linguistic_complexity(string):
    substring_count = get_substring_count(string)
    maximum_substring_count = get_maximum_substring_count(string)
    return float(substring_count) / float(maximum_substring_count)
        
def get_longest_common_prefix(x, y):
    n = min(len(x), len(y))
    for i in range(n):
        if x[i] != y[i]:
            return i
    return n
    
def main():
    dna_string = sys.stdin.readline().strip()
    print(get_linguistic_complexity(dna_string))
    
if __name__ == '__main__':
    sys.exit(main())
