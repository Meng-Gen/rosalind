import math
import sys

def C(m, n):
    return math.factorial(m) // math.factorial(n) // math.factorial(m-n)
    
def main():
    n = int(sys.stdin.readline().strip())
    array_size = 2*n
    prob_list = [C(array_size, i) * 0.5**(array_size) for i in range(array_size + 1)]
    acc_prob_list = [prob_list[0]]
    for i in range(1, len(prob_list)):
        acc_prob_list.append(acc_prob_list[-1] + prob_list[i])
    prob_common_log_list = [math.log10(_) for _ in acc_prob_list[::-1]][1:]
    print(' '.join(map(str, prob_common_log_list)))
    
if __name__ == '__main__':
    sys.exit(main())
