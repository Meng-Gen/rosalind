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

class Blosum62ScoringMatrix():
    def __init__(self):
        self.blosum62_column = {}
        self.blosum62_table = [   
            [ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2],
            [ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
            [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3],
            [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2],
            [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3],
            [ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3],
            [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2],
            [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1],
            [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2],
            [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1],
            [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1],
            [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2],
            [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3],
            [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1],
            [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2],
            [ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2],
            [ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2],
            [ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1],
            [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2],
            [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7],
        ]
        self.__init_blosum62_column()

    def get_score(self, x, y):
        return self.blosum62_table[self.blosum62_column[x]][self.blosum62_column[y]]
    
    def __init_blosum62_column(self):
        blosum62_column_name = [
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
        ] 
        for i in range(len(blosum62_column_name)):
            self.blosum62_column[blosum62_column_name[i]] = i
    
class GlobalAlignment():
    def __init__(self, gap_opening_penalty, gap_extension_penalty):
        self.score_matrix = Blosum62ScoringMatrix()
        self.h = gap_opening_penalty - gap_extension_penalty
        self.g = gap_extension_penalty
        
    def edit_distance(self, s, t):
        m, n = len(s), len(t)
        M = [[-sys.maxsize for _ in range(n+1)] for _ in range(m+1)] 
        I = [[-sys.maxsize for _ in range(n+1)] for _ in range(m+1)] 
        J = [[-sys.maxsize for _ in range(n+1)] for _ in range(m+1)] 
        M[0][0] = 0
        for i in range(m+1):
            I[i][0] = self.h + self.g * i
        for j in range(n+1):
            J[0][j] = self.h + self.g * j
        for j in range(n):
            for i in range(m):
                score = self.score_matrix.get_score(s[i], t[j])
                M[i+1][j+1] = max(M[i][j] + score, I[i][j] + score, J[i][j] + score)
                I[i+1][j+1] = max(M[i][j+1] + self.h + self.g, I[i][j+1] + self.g)
                J[i+1][j+1] = max(M[i+1][j] + self.h + self.g, J[i+1][j] + self.g)
        best_score = max(M[m][n], I[m][n], J[m][n])
    
        # Traceback
        i, j = m, n
        s_traceback, t_traceback = '', ''
        traceback_table = None
        if M[m][n] == best_score:
            traceback_table = [M]
        elif I[m][n] == best_score:
            traceback_table = [I]
        elif J[m][n] == best_score:
            traceback_table = [J]
        
        while i != 0 or j != 0:
            if M in traceback_table:
                s_traceback += s[i-1]
                t_traceback += t[j-1]
                score = self.score_matrix.get_score(s[i-1], t[j-1])
                if M[i-1][j-1] + score == M[i][j]:
                    traceback_table = [M]
                elif I[i-1][j-1] + score == M[i][j]:
                    traceback_table = [I]
                elif J[i-1][j-1] + score == M[i][j]:
                    traceback_table = [J]
                i, j = i-1, j-1
            elif I in traceback_table:
                s_traceback += s[i-1]
                t_traceback += '-'
                if M[i-1][j] + self.h + self.g == I[i][j]:
                    traceback_table = [M]
                elif I[i-1][j] + self.g == I[i][j]:
                    traceback_table = [I]
                i, j = i-1, j
            elif J in traceback_table:
                s_traceback += '-'
                t_traceback += t[j-1]
                if M[i][j-1] + self.h + self.g == J[i][j]:
                    traceback_table = [M]
                elif J[i][j-1] + self.g == J[i][j]:
                    traceback_table = [J]
                i, j = i, j-1
        return best_score, s_traceback[::-1], t_traceback[::-1]
    
def main():
    fasta_records = read_fasta_records()
    s, t = fasta_records[0][1], fasta_records[1][1]
    alignment = GlobalAlignment(-11, -1)
    print('\n'.join(map(str, alignment.edit_distance(s, t))))

if __name__ == '__main__':
    sys.exit(main())
