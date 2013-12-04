import math
import sys

def C(m, n):
    return math.factorial(m) // math.factorial(n) // math.factorial(m-n)    

def matrix_multiply(A, B):
    A_row_dim, A_col_dim = len(A), len(A[0])
    B_row_dim, B_col_dim = len(B), len(B[0])
    assert(A_col_dim == B_row_dim)
    C = [[0.0 for row in range(B_col_dim)] for col in range(A_row_dim)]
    for i in range(A_row_dim):
        for j in range(B_col_dim):
            for k in range(A_col_dim):
                C[i][j] += A[i][k] * B[k][j]
    return C
    
class WrightFisherModel():
    def __init__(self, N, m):
        self.N = N
        self.p = m/(2*N)
        self.q = 1 - m/(2*N)
        self.theta = None
        self.T = [[None for i in range(2*N+1)] for j in range(2*N+1)]

    def build(self, generation):
        self.__build_theta(self.N, self.p, self.q)
        self.__build_T(self.N)
        curr_theta = self.theta
        for g in range(generation - 1):
            curr_theta = matrix_multiply(self.T, curr_theta)
        return curr_theta
        
    def __build_theta(self, N, p, q):
        self.theta = [[C(2*N, j) * p**j * q**(2*N-j)] for j in range(2*N + 1)]
        
    def __build_T(self, N):
        for i in range(2*N + 1):
            for j in range(2*N + 1):
                self.T[j][i] = C(2*N, j) * (i/(2*N))**j * (1-i/(2*N))**(2*N-j)
    
def get_final_prob(theta, k):
    theta_dim = len(theta)
    sum = 0.0    
    for i in range(theta_dim - k):    
        sum += theta[i][0]
    return sum
    
def main():
    N, m, g, k = map(int, sys.stdin.readline().strip().split())
    model = WrightFisherModel(N, m)
    theta = model.build(g)
    print(get_final_prob(theta, k))
        
if __name__ == '__main__':
    sys.exit(main())
