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
    
class OverlapAlignment():
    def __init__(self):
        self.g = -1
        
    def edit_distance(self, s, t):
        m, n = len(s), len(t)
        M = [[0 for _ in range(n+1)] for _ in range(m+1)] 
        for i in range(m):
            M[i+1][0] = self.g * (i+1)
        for j in range(n):
            M[0][j+1] = self.g * (j+1)
        for j in range(n):
            for i in range(m):
                score = self.get_score(s[i], t[j])
                M[i+1][j+1] = max(M[i+1][j] + self.g, M[i][j+1] + self.g, M[i][j] + score)
        # traceback
        s_traceback, t_traceback = '', ''
        i, j = m, n
        while i != 0 and j != 0:
            if M[i][j-1] + self.g == M[i][j]:
                s_traceback += '-'
                t_traceback += t[j-1]
                i, j = i, j-1
            elif M[i-1][j] + self.g == M[i][j]:
                s_traceback += s[i-1]
                t_traceback += '-'
                i, j = i-1, j
            else:
                s_traceback += s[i-1]
                t_traceback += t[j-1]
                i, j = i-1, j-1
        return M[m][n], s_traceback[::-1], t_traceback[::-1]
    
    def get_score(self, x, y):
        if x == y:
            return 1
        else:
            return -30000
    
def main():
    fasta_records = read_fasta_records()
    s, t = fasta_records[0][1], fasta_records[1][1]
    alignment = OverlapAlignment()
    score, s_alignment, t_alignment = alignment.edit_distance(s, t)
    print(s_alignment.count('-') + t_alignment.count('-'))

if __name__ == '__main__':
    sys.exit(main())
