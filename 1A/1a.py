import sys

dna = sys.stdin.readline().strip()
k = int(sys.stdin.readline().strip())

freq = {}
for i in range(len(dna) + 1 - k):
    mer = dna[i:i+k]
    if not mer in freq:
        freq[mer] = 0
    freq[mer] += 1

max_freq = 0
max_freq_mer = []
for mer in freq:
    if freq[mer] > max_freq:
        max_freq = freq[mer]
        max_freq_mer = []
    if freq[mer] == max_freq:
        max_freq_mer.append(mer)
print(' '.join(max_freq_mer))
