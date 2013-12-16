import sys

def get_character_table(dna_strings):
    character_table = []
    m = len(dna_strings[0])
    for j in range(m):
        character = get_character(dna_strings, j)
        if character:
            character_table.append(character)
    return character_table

def get_character(dna_strings, j):
    n = len(dna_strings)
    used_symbol = set()
    dict = { 'A': set(), 'C': set(), 'T': set(), 'G': set() }
    for i in range(n):
        symbol = dna_strings[i][j]
        used_symbol.add(symbol)
        dict[symbol].add(i)

    # Cannot create character
    if len(used_symbol) != 2:
        return None

    x, y = used_symbol
    # Ignore trivial character
    if len(dict[x]) == 1 or len(dict[y]) == 1:
        return None

    character = [1 for _ in range(n)]
    for i in dict[x]:
        character[i] = 0
    return ''.join(map(str, character))
        
def main():
    dna_strings = [line.strip() for line in sys.stdin.readlines()]
    print('\n'.join(get_character_table(dna_strings)))
    
if __name__ == '__main__':
    sys.exit(main())
