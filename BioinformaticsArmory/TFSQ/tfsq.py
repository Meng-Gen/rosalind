import sys

def main():
    from Bio import SeqIO
    fastq_dict = SeqIO.index("in.txt", "fastq")
    for fastq in fastq_dict:
        print(fastq_dict[fastq].format("fasta"))    
    
if __name__ == '__main__':
    sys.exit(main())
