import sys

from Bio import ExPASy
from Bio import SwissProt

def get_records(ids):
    records = []
    for id in ids:
        handle = ExPASy.get_sprot_raw(id)
        record = SwissProt.read(handle)
        records.append(record.sequence)
    return records
    
def main():
    ids = sys.stdin.read().strip().split()
    records = get_records(ids)
    print('\n'.join(records))
    
    # Then go to http://www.ebi.ac.uk/Tools/psa/emboss_water/
    
if __name__ == '__main__':
    sys.exit(main())
