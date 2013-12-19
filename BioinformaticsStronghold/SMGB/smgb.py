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
    
class SemiglobalAlignment():
    def __init__(self):
        self.g = -1
        
    def edit_distance(self, s, t):
        m, n = len(s), len(t)
        M = [[0 for _ in range(n+1)] for _ in range(m+1)] 
        for j in range(n):
            for i in range(m):
                score = self.get_score(s[i], t[j])
                M[i+1][j+1] = max(M[i+1][j] + self.g, M[i][j+1] + self.g, M[i][j] + score)
        # get best score
        best_score_so_far = M[0][n]
        best_pos_so_far = None
        for i in range(m):
            if M[i+1][n] > best_score_so_far:
                best_score_so_far = M[i+1][n]
                best_pos_so_far = i+1, n
        for j in range(n):
            if M[m][j+1] > best_score_so_far:
                best_score_so_far = M[m][j+1]
                best_pos_so_far = m, j+1
        # traceback
        s_traceback = s[best_pos_so_far[0]:][::-1] + '-' * (n - best_pos_so_far[1])
        t_traceback = t[best_pos_so_far[1]:][::-1] + '-' * (m - best_pos_so_far[0])
        i, j = best_pos_so_far
        while i != 0 or j != 0:
            if i == 0:
                s_traceback += ('-' * j)
                t_traceback += t[:j][::-1]
                break
            elif j == 0:
                s_traceback += s[:i][::-1]
                t_traceback += ('-' * i)
                break
            elif M[i][j-1] + self.g == M[i][j]:
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
        return best_score_so_far, s_traceback[::-1], t_traceback[::-1]
    
    def get_score(self, x, y):
        if x == y:
            return 1
        else:
            return -1
    
def main():
    fasta_records = read_fasta_records()
    s, t = fasta_records[0][1], fasta_records[1][1]
    alignment = SemiglobalAlignment()
    score, s_alignment, t_alignment = alignment.edit_distance(s, t)
    print(score)
    print(s_alignment)
    print(t_alignment)

if __name__ == '__main__':
    sys.exit(main())
