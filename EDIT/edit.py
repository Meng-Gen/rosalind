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

def levenshtein_distance(s, t):
    m, n = len(s), len(t)
    d = [[0 for _ in range(n+1)] for _ in range(m+1)] 
    
    for i in range(m):
        d[i+1][0] = i+1
    for j in range(n):
        d[0][j+1] = j+1

    for j in range(n):
        for i in range(m):
            if s[i] == t[j]:
                d[i+1][j+1] = d[i][j]
            else:
                d[i+1][j+1] = min(d[i+1][j], d[i][j+1], d[i][j]) + 1
    return d[m][n]

def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    s, t = fasta_records[0][1], fasta_records[1][1]
    print(levenshtein_distance(s, t))

if __name__ == '__main__':
    sys.exit(main())
