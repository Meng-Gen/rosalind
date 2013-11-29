import sys
import urllib.request

def get_protein_data(uniprot_id):
    url = '''http://www.uniprot.org/uniprot/{uniprot_id}.txt'''.format(uniprot_id=uniprot_id)
    response = urllib.request.urlopen(url)
    return response.read().decode('utf-8').split('\n')
        
def find_biological_processes(uniprot_id):
    protein_data = get_protein_data(uniprot_id)
    for line in protein_data:
        # find gene ontology (GO) section
        if not line.startswith('DR'):
            continue
        details = [_.strip() for _ in line[2:].split(';')]
        # find gene ontology (GO) section
        if details[0] != 'GO':
            continue
        # find biological processes
        if not details[2].startswith('P:'):
            continue
        print(details[2][2:])
        
def main():
    for uniprot_id in sys.stdin.readlines():
        find_biological_processes(uniprot_id.strip())
    
if __name__ == '__main__':
    sys.exit(main())
