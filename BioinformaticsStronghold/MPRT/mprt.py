import sys
import urllib.request

def get_fasta_records(lines):
    fasta_records = []
    is_first_record = True
    id, record = '', ''
    for line in lines:
        line = line.strip();
        if line and line[0] == '>':
            if is_first_record:
                is_first_record = False
            else:
                fasta_records.append([id, record])
            id, record = line[1:], ''
        else:
            record += line
    fasta_records.append([id, record])
    return fasta_records

def get_protein_sequence(uniprot_id):
    url = '''http://www.uniprot.org/uniprot/{uniprot_id}.fasta'''.format(uniprot_id=uniprot_id)
    response = urllib.request.urlopen(url)
    return response.read().decode('utf-8')

def n_glycosylation_positions(protein):
    positions = [] # 1-based position
    for i in range(len(protein)+1-4):
        motif = protein[i:i+4]
        if motif[0] == 'N' and motif[1] != 'P' and motif[2] in ['S', 'T'] and motif[3] != 'P':
            positions.append(i+1) 
    return positions

def find_protein_motif(uniprot_id):
    lines = get_protein_sequence(uniprot_id).split('\n')
    record = get_fasta_records(lines)[0]
    positions = n_glycosylation_positions(record[1])
    if positions:
        print(uniprot_id)
        print(' '.join(map(str, positions)))
        
def main():
    for uniprot_id in sys.stdin.readlines():
        find_protein_motif(uniprot_id.strip())
    
if __name__ == '__main__':
    sys.exit(main())
