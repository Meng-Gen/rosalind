import sys

def get_fasta_records(lines):
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

def edit_distance(s, t):
    m, n = len(s), len(t)
    d = [[0 for _ in range(n+1)] for _ in range(m+1)] 
    a = [['?' for _ in range(n+1)] for _ in range(m+1)] 
    
    for i in range(m):
        d[i+1][0] = i+1
    for j in range(n):
        d[0][j+1] = j+1
    a[0][0] = '+'
        
    for j in range(n):
        for i in range(m):
            if s[i] == t[j]:
                d[i+1][j+1] = d[i][j]
                a[i+1][j+1] = '+'
            else:
                d[i+1][j+1] = min(d[i+1][j], d[i][j+1], d[i][j]) + 1
                if d[i+1][j+1] == d[i][j] + 1:
                    a[i+1][j+1] = '#'
                if d[i+1][j+1] == d[i][j+1] + 1:
                    a[i+1][j+1] = '|'
                if d[i+1][j+1] == d[i+1][j] + 1:
                    a[i+1][j+1] = '-'
    
    s_reversed_alignment, t_reversed_alignment = '', ''
    i, j = m, n
    while i != 0 and j != 0:
        if a[i][j] == '-':
            s_reversed_alignment += '-'
            t_reversed_alignment += t[j-1]
            i, j = i, j-1
        elif a[i][j] == '|':
            s_reversed_alignment += s[i-1]
            t_reversed_alignment += '-'
            i, j = i-1, j
        elif a[i][j] == '#' or a[i][j] == '+':
            s_reversed_alignment += s[i-1]
            t_reversed_alignment += t[j-1]
            i, j = i-1, j-1
    return d[m][n], s_reversed_alignment[::-1], t_reversed_alignment[::-1]
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    s, t = fasta_records[0][1], fasta_records[1][1]
    print('\n'.join(map(str, edit_distance(s, t))))

if __name__ == '__main__':
    sys.exit(main())
