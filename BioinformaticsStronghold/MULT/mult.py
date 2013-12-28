import itertools
import sys

def read_fasta_records():
    lines = sys.stdin.readlines()
    records = []
    is_first_record = True
    id, content = '', ''
    for line in lines:
        line = line.strip();
        if line[0] == '>':
            if is_first_record:
                is_first_record = False
            else:
                records.append([id, content])
            id, content = line[1:], ''
        else:
            content += line
    records.append([id, content])
    return records

class GlobalAlignment():
    def __init__(self, sequence_array):
        self.sequence_array = sequence_array
        self.M = None
               
    def edit_distance(self):
        self.init_M()
        self.set_0_dim()
        self.set_1_dim()
        self.set_2_dim()
        self.set_3_dim()
        self.set_4_dim()
        return self.traceback()
    
    def init_M(self):
        self.M = {}
        range_array = [range(len(_)+1) for _ in self.sequence_array]
        for index_set in itertools.product(*range_array):
            self.M[self.hash(index_set)] = -1000

    def set_0_dim(self):
        n = len(self.sequence_array)
        self.M[self.hash([0 for _ in range(n)])] = 0

    def set_1_dim(self):
        range_array = [range(len(_)+1) for _ in self.sequence_array]
        n = len(range_array)
        for i in range(n):
            index_set = [0 for _ in range(n)]
            for j in range_array[i]:
                index_set[i] = j
                self.M[self.hash(index_set)] = -3 * j
    
    def set_2_dim(self):
        range_array = [range(1, len(_)+1) for _ in self.sequence_array]
        n = len(range_array)
        for i_pos, j_pos in itertools.combinations(range(n), 2):
            index_set = [0 for _ in range(n)]
            index_i_set = [0 for _ in range(n)]
            index_j_set = [0 for _ in range(n)]
            index_ij_set = [0 for _ in range(n)]
            s, t = self.sequence_array[i_pos], self.sequence_array[j_pos]
            for i, j in itertools.product(range_array[i_pos], range_array[j_pos]):
                index_set[i_pos], index_set[j_pos] = i, j
                index_i_set[i_pos], index_i_set[j_pos] = i-1, j
                index_j_set[i_pos], index_j_set[j_pos] = i, j-1
                index_ij_set[i_pos], index_ij_set[j_pos] = i-1, j-1
                score = self.get_score(s[i-1], t[j-1])
                
                self.M[self.hash(index_set)] = max( \
                    self.get_M(index_i_set) - 1, \
                    self.get_M(index_j_set) - 1, \
                    self.get_M(index_ij_set) + score, \
                )
            for i, j in itertools.product(range_array[i_pos], range_array[j_pos]):
                index_set[i_pos], index_set[j_pos] = i, j
                self.M[self.hash(index_set)] -= 2 * (i+j)
    
    def set_3_dim(self):
        range_array = [range(1, len(_)+1) for _ in self.sequence_array]
        n = len(range_array)
        for i_pos, j_pos, k_pos in itertools.combinations(range(n), 3):
            index_set = [0 for _ in range(n)]
            index_i_set = [0 for _ in range(n)]
            index_j_set = [0 for _ in range(n)]
            index_k_set = [0 for _ in range(n)]
            index_ij_set = [0 for _ in range(n)]
            index_ik_set = [0 for _ in range(n)]
            index_jk_set = [0 for _ in range(n)]
            index_ijk_set = [0 for _ in range(n)]
            s = self.sequence_array[i_pos]
            t = self.sequence_array[j_pos]
            u = self.sequence_array[k_pos]
            for i, j, k in itertools.product(range_array[i_pos], range_array[j_pos], range_array[k_pos]):
                index_set[i_pos], index_set[j_pos], index_set[k_pos] = i, j, k
                index_i_set[i_pos], index_i_set[j_pos], index_i_set[k_pos] = i-1, j, k
                index_j_set[i_pos], index_j_set[j_pos], index_j_set[k_pos] = i, j-1, k
                index_k_set[i_pos], index_k_set[j_pos], index_k_set[k_pos] = i, j, k-1
                index_ij_set[i_pos], index_ij_set[j_pos], index_ij_set[k_pos] = i-1, j-1, k
                index_ik_set[i_pos], index_ik_set[j_pos], index_ik_set[k_pos] = i-1, j, k-1
                index_jk_set[i_pos], index_jk_set[j_pos], index_jk_set[k_pos] = i, j-1, k-1
                index_ijk_set[i_pos], index_ijk_set[j_pos], index_ijk_set[k_pos] = i-1, j-1, k-1
                score_ij = self.get_score(s[i-1], t[j-1])
                score_ik = self.get_score(s[i-1], u[k-1])
                score_jk = self.get_score(t[j-1], u[k-1])
                
                self.M[self.hash(index_set)] = max( \
                    self.get_M(index_i_set) - 2, \
                    self.get_M(index_j_set) - 2, \
                    self.get_M(index_k_set) - 2, \
                    self.get_M(index_ij_set) + score_ij - 2, \
                    self.get_M(index_ik_set) + score_ik - 2, \
                    self.get_M(index_jk_set) + score_jk - 2, \
                    self.get_M(index_ijk_set) + score_ij + score_ik + score_jk, \
                )
            for i, j, k in itertools.product(range_array[i_pos], range_array[j_pos], range_array[k_pos]):
                index_set[i_pos], index_set[j_pos], index_set[k_pos] = i, j, k
                self.M[self.hash(index_set)] -= (i+j+k)
    
    def set_4_dim(self):
        range_array = [range(1, len(_)+1) for _ in self.sequence_array]
        s, t, u, v = self.sequence_array
        for i, j, k, l in itertools.product(*range_array):
            score_ij = self.get_score(s[i-1], t[j-1])
            score_ik = self.get_score(s[i-1], u[k-1])
            score_il = self.get_score(s[i-1], v[l-1])
            score_jk = self.get_score(t[j-1], u[k-1])
            score_jl = self.get_score(t[j-1], v[l-1])
            score_kl = self.get_score(u[k-1], v[l-1])
        
            self.M[self.hash([i, j, k, l])] = max( \
                self.get_M([i-1, j, k, l]) - 3, \
                self.get_M([i, j-1, k, l]) - 3, \
                self.get_M([i, j, k-1, l]) - 3, \
                self.get_M([i, j, k, l-1]) - 3, \
                self.get_M([i-1, j-1, k, l]) + score_ij - 4, \
                self.get_M([i-1, j, k-1, l]) + score_ik - 4, \
                self.get_M([i-1, j, k, l-1]) + score_il - 4, \
                self.get_M([i, j-1, k-1, l]) + score_jk - 4, \
                self.get_M([i, j-1, k, l-1]) + score_jl - 4, \
                self.get_M([i, j, k-1, l-1]) + score_kl - 4, \
                self.get_M([i-1, j-1, k-1, l]) + score_ij + score_ik + score_jk - 3, \
                self.get_M([i-1, j-1, k, l-1]) + score_ij + score_il + score_jl - 3, \
                self.get_M([i-1, j, k-1, l-1]) + score_ik + score_il + score_kl - 3, \
                self.get_M([i, j-1, k-1, l-1]) + score_jk + score_jl + score_kl - 3, \
                self.get_M([i-1, j-1, k-1, l-1]) + score_ij + score_ik + score_il + \
                                                   score_jk + score_jl + score_kl, \
            )     
    
    def traceback(self):
        assert(len(self.sequence_array) == 4)
        s, t, u, v = self.sequence_array
        i, j, k, l = [len(_) for _ in self.sequence_array]
        best_score = self.get_M([i, j, k, l])
        s_traceback, t_traceback, u_traceback, v_traceback = '', '', '', ''
        
        while i != 0 or j != 0 or k != 0 or l != 0:
            score = self.get_M([i, j, k, l])
            score_ij = self.get_score(s[i-1], t[j-1])
            score_ik = self.get_score(s[i-1], u[k-1])
            score_il = self.get_score(s[i-1], v[l-1])
            score_jk = self.get_score(t[j-1], u[k-1])
            score_jl = self.get_score(t[j-1], v[l-1])
            score_kl = self.get_score(u[k-1], v[l-1])
            if self.get_M([i-1, j, k, l]) - 3 == score:
                s_traceback += s[i-1]
                t_traceback += '-'
                u_traceback += '-'
                v_traceback += '-'
                i, j, k, l = i-1, j, k, l
            elif self.get_M([i, j-1, k, l]) - 3 == score:
                s_traceback += '-'
                t_traceback += t[j-1]
                u_traceback += '-'
                v_traceback += '-'
                i, j, k, l = i, j-1, k, l
            elif self.get_M([i, j, k-1, l]) - 3 == score:
                s_traceback += '-'
                t_traceback += '-'
                u_traceback += u[k-1]
                v_traceback += '-'
                i, j, k, l = i, j, k-1, l
            elif self.get_M([i, j, k, l-1]) - 3 == score:
                s_traceback += '-'
                t_traceback += '-'
                u_traceback += '-'
                v_traceback += v[l-1]
                i, j, k, l = i, j, k, l-1
            elif self.get_M([i-1, j-1, k, l]) + score_ij - 4 == score:
                s_traceback += s[i-1]
                t_traceback += t[j-1]
                u_traceback += '-'
                v_traceback += '-'
                i, j, k, l = i-1, j-1, k, l
            elif self.get_M([i-1, j, k-1, l]) + score_ik - 4 == score:
                s_traceback += s[i-1]
                t_traceback += '-'
                u_traceback += u[k-1]
                v_traceback += '-'
                i, j, k, l = i-1, j, k-1, l            
            elif self.get_M([i-1, j, k, l-1]) + score_il - 4 == score:
                s_traceback += s[i-1]
                t_traceback += '-'
                u_traceback += '-'
                v_traceback += v[l-1]
                i, j, k, l = i-1, j, k, l-1   
            elif self.get_M([i, j-1, k-1, l]) + score_jk - 4 == score:
                s_traceback += '-'
                t_traceback += t[j-1]
                u_traceback += u[k-1]
                v_traceback += '-'
                i, j, k, l = i, j-1, k-1, l
            elif self.get_M([i, j-1, k, l-1]) + score_jl - 4 == score:
                s_traceback += '-'
                t_traceback += t[j-1]
                u_traceback += '-'
                v_traceback += v[l-1]
                i, j, k, l = i, j-1, k, l-1
            elif self.get_M([i, j, k-1, l-1]) + score_kl - 4 == score:
                s_traceback += '-'
                t_traceback += '-'
                u_traceback += u[k-1]
                v_traceback += v[l-1]
                i, j, k, l = i, j, k-1, l-1
            elif self.get_M([i-1, j-1, k-1, l]) + score_ij + score_ik + score_jk - 3 == score:
                s_traceback += s[i-1]
                t_traceback += t[j-1]
                u_traceback += u[k-1]
                v_traceback += '-'            
                i, j, k, l = i-1, j-1, k-1, l
            elif self.get_M([i-1, j-1, k, l-1]) + score_ij + score_il + score_jl - 3 == score:
                s_traceback += s[i-1]
                t_traceback += t[j-1]
                u_traceback += '-'
                v_traceback += v[l-1]            
                i, j, k, l = i-1, j-1, k, l-1
            elif self.get_M([i-1, j, k-1, l-1]) + score_ik + score_il + score_kl - 3 == score:
                s_traceback += s[i-1]
                t_traceback += '-'
                u_traceback += u[k-1]
                v_traceback += v[l-1]            
                i, j, k, l = i-1, j, k-1, l-1
            elif self.get_M([i, j-1, k-1, l-1]) + score_jk + score_jl + score_kl - 3 == score:
                s_traceback += '-'
                t_traceback += t[j-1]
                u_traceback += u[k-1]
                v_traceback += v[l-1]
                i, j, k, l = i, j-1, k-1, l-1
            else:
                s_traceback += s[i-1]
                t_traceback += t[j-1]
                u_traceback += u[k-1]
                v_traceback += v[l-1]
                i, j, k, l = i-1, j-1, k-1, l-1
        return best_score, s_traceback[::-1], t_traceback[::-1], u_traceback[::-1], v_traceback[::-1]
    
    def get_M(self, index):
        h = self.hash(index)
        if h in self.M:
            return self.M[h]
        else:
            return -1000
    
    def debug_M(self):
        range_array = [range(len(_)+1) for _ in self.sequence_array]
        for index_set in itertools.product(*range_array):
            print(index_set, '=>', self.M[self.hash(index_set)])
    
    def get_score(self, x, y):
        if x == y:
            return 0
        else:
            return -1
    
    def hash(self, index):
        return '|'.join(map(str, index))
    
def main():
    fasta_records = read_fasta_records()
    sequence_array = [fasta_records[i][1] for i in range(4)]
    alignment = GlobalAlignment(sequence_array)
    print('\n'.join(map(str, alignment.edit_distance())))

if __name__ == '__main__':
    sys.exit(main())
