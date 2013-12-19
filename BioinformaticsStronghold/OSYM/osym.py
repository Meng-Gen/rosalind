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
    def __init__(self):
        self.g = -1
        
    def edit_distance(self, s, t):
        m, n = len(s), len(t)
        A = self.get_score_matrix(s, t)
        B = self.get_score_matrix(s[::-1], t[::-1])
        M = [[None for j in range(n)] for i in range(m)] 
        for j in range(n):
            for i in range(m):
                M[i][j] = A[i][j] + self.get_score(s[i], t[j]) + B[m-1-i][n-1-j]
        total = sum(sum(_) for _ in M)
        return A[m][n], total
    
    def get_score_matrix(self, s, t):
        m, n = len(s), len(t)
        d = [[0 for _ in range(n+1)] for _ in range(m+1)] 
        for i in range(m):
            d[i+1][0] = self.g * (i+1)
        for j in range(n):
            d[0][j+1] = self.g * (j+1)
        for j in range(n):
            for i in range(m):
                score = self.get_score(s[i], t[j])    
                d[i+1][j+1] = max(d[i+1][j] + self.g, d[i][j+1] + self.g, d[i][j] + score)
        return d
        
    def get_score(self, x, y):
        if x == y:
            return 1
        else:
            return -1
    
def main():
    fasta_records = read_fasta_records()
    s, t = fasta_records[0][1], fasta_records[1][1]
    alignment = GlobalAlignment()
    print('\n'.join(map(str, alignment.edit_distance(s, t))))

if __name__ == '__main__':
    sys.exit(main())
