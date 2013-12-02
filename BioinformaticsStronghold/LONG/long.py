import sys

def get_fasta_records(lines):
    dna_strings = []
    is_first_dna = True
    id, dna = '', ''
    for line in lines:
        line = line.strip();
        if line[0] == '>':
            if is_first_dna:
                is_first_dna = False
            else:
                dna_strings.append([id, dna])
            id, dna = line[1:], ''
        else:
            dna += line
    dna_strings.append([id, dna])
    return dna_strings

def get_shortest_superstring(fasta_records):
    adjacent_matrix = get_adjacent_matrix(fasta_records)
    assembly_chain = get_assembly_chain(adjacent_matrix)
    return glue(fasta_records, assembly_chain)
    
def get_adjacent_matrix(fasta_records):
    num_vertex = len(fasta_records)
    return [[is_adjacent(fasta_records[x][1], fasta_records[y][1]) 
        for y in range(num_vertex)] for x in range(num_vertex)]

def is_adjacent(u, v):
    for i in range(len(u)//2):
        j = len(u) - i
        if u[i:] == v[:j] and j > len(v)//2:
            return True
    return False        
        
def get_assembly_chain(adjacent_matrix):
    num_vertex = len(adjacent_matrix)
    reversed_chain = [get_final_vertex(adjacent_matrix)]
    while len(reversed_chain) != num_vertex:
        last_vertex = reversed_chain[-1]
        for i in range(num_vertex):
            if adjacent_matrix[i][last_vertex] is True and i != last_vertex:
                reversed_chain.append(i)
                break
    return reversed_chain[::-1]

def get_final_vertex(adjacent_matrix):
    num_vertex = len(adjacent_matrix)
    final_vertex_list = [i for i in range(num_vertex) 
        if adjacent_matrix[i].count(False) == num_vertex - 1]
    assert(len(final_vertex_list) == 1)
    return final_vertex_list[0]    
    
def glue(fasta_records, assembly_chain):
    superstring = fasta_records[assembly_chain[0]][1]
    num_vertex = len(fasta_records)
    for i in range(1, num_vertex):
        prev_vertex = fasta_records[assembly_chain[i-1]][1]
        curr_vertex = fasta_records[assembly_chain[i]][1]
        for j in range(len(curr_vertex), 0, -1):
            if curr_vertex[:j] == prev_vertex[len(prev_vertex) - j:]:
                superstring += curr_vertex[j:]
                break
    return superstring
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    print(get_shortest_superstring(fasta_records))

if __name__ == '__main__':
    sys.exit(main())
