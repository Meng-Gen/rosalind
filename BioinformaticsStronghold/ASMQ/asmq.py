import math
import sys
      
def get_max_L(dna_length_array, n_statistic):
    dna_length_array.sort()
    dna_length_sum = sum(dna_length_array)
    threshold = math.ceil(dna_length_sum * n_statistic)    
    curr_sum = 0
    for length in dna_length_array[::-1]:
        curr_sum += length
        if curr_sum >= threshold:
            return length
      
def main():
    dna_length_array = [len(dna.strip()) for dna in sys.stdin.readlines()]
    print(get_max_L(dna_length_array, 0.50), get_max_L(dna_length_array, 0.75))
    
if __name__ == '__main__':
    sys.exit(main())
