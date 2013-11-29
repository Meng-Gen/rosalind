import sys

freq = {}
for word in sys.stdin.read().strip().split():
    if word not in freq:
        freq[word] = 0
    freq[word] += 1
for word in freq:
    print(word, freq[word])
