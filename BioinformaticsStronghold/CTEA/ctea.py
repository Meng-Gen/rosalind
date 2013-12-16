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

def edit_distance_ways(s, t):
    m, n = len(s), len(t)
    d = [[0 for _ in range(n+1)] for _ in range(m+1)] 
    ways = [[0 for _ in range(n+1)] for _ in range(m+1)] 
    for i in range(m+1):
        d[i][0] = i
        ways[i][0] = 1
    for j in range(n+1):
        d[0][j] = j
        ways[0][j] = 1
        
    for j in range(1, n+1):
        for i in range(1, m+1):
            if s[i-1] == t[j-1]:
                d[i][j] = d[i-1][j-1]
                ways[i][j] += ways[i-1][j-1]
            else:                
                d[i][j] = min(d[i-1][j-1] + 1, d[i-1][j] + 1, d[i][j-1] + 1)
            if d[i][j] == d[i-1][j-1] + 1:
                ways[i][j] += ways[i-1][j-1]
            if d[i][j] == d[i-1][j] + 1:
                ways[i][j] += ways[i-1][j]        
            if d[i][j] == d[i][j-1] + 1:
                ways[i][j] += ways[i][j-1]

    return ways[m][n]
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    s, t = fasta_records[0][1], fasta_records[1][1]
    print(edit_distance_ways(s, t) % 134217727)

if __name__ == '__main__':
    sys.exit(main())
