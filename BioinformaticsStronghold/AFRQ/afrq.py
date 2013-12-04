import math
import sys

def get_disease_carriers_prob(x):
    q = math.sqrt(x)
    p = 1 - q
    return 2*p*q + q*q

def main():
    homozygous_recessive_prob_list = list(map(float, sys.stdin.readline().strip().split()))
    rv = [get_disease_carriers_prob(x) for x in homozygous_recessive_prob_list]
    print(' '.join(map(str, rv)))
    
if __name__ == '__main__':
    sys.exit(main())
