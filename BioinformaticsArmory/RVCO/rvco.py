import sys

def main():
    from Bio import SeqIO
    count = 0
    for record in SeqIO.parse('in.txt', 'fasta'):
        if record.seq.tostring() == record.reverse_complement().seq.tostring():
            count += 1
    print(count)
    
if __name__ == '__main__':
    sys.exit(main())
