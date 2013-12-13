import sys

from Bio import SeqIO

def read_quality_threshold():
    with open('in.txt') as dataset_file:
        line = dataset_file.readline()
        return int(line.strip())
        
def read_records():
    return SeqIO.parse('in.txt', 'fastq')

def get_average_quality(record):
    phred_quality_array = record.letter_annotations["phred_quality"]
    return float(sum(phred_quality_array)) / float(len(phred_quality_array))
    
def main():
    quality_threshold = read_quality_threshold()
    count = 0
    for record in read_records():
        if get_average_quality(record) < quality_threshold:
            count += 1
    print(count)
    
if __name__ == '__main__':
    sys.exit(main())
