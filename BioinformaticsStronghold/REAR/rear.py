import ast
import itertools
import os.path
import sys

def read_permutation_pairs():
    pairs = []
    is_first_permutation = True
    p, q = None, None
    for line in sys.stdin.readlines():
        line = line.strip()
        if not line:
            continue
        if is_first_permutation:
            p = list(map(int, line.split()))
        else:
            q = list(map(int, line.split()))
            pairs.append([p, q])
        is_first_permutation = not is_first_permutation
    return pairs

def normalize(pair):
    p, q, r = pair[0], pair[1], []
    for i in range(len(q)):
        for j in range(len(p)):
            if p[j] == q[i]:
                r.append(j+1)
                break
    return r

class ReversalDistanceQuery():
    def __init__(self, n):
        self.n = n
        self.distance_map_data_file = 'DistanceMapData[%d].txt'
        self.distance_map = {}
        self.__init_distance_map()
        self.__precalculate()
        
    def __write_data_file(self, distance):
        file = self.distance_map_data_file % distance
        with open(file, 'w') as file_handle:
            file_handle.write(str(self.distance_map[distance]))
        
    def __read_data_file(self, distance):
        file = self.distance_map_data_file % distance
        with open(file, 'r') as file_handle:
            self.distance_map[distance] = ast.literal_eval(file_handle.read())
        
    def __has_data_file(self, distance):
        file = self.distance_map_data_file % distance
        return os.path.isfile(file)
        
    def __precalculate(self):
        for i in range(2, self.n):
            self.__update(i)
            print(len(self.distance_map[i]))
    
    def __init_distance_map(self):
        for i in range(self.n):
            self.distance_map[i] = set()
        initial_permutation = [_ for _ in range(1, self.n+1)]
        self.__insert(initial_permutation, 0)
        for reversal in self.__reversal_generator(initial_permutation):
            self.__insert(reversal, 1)
    
    def __update(self, distance):
        if self.__has_data_file(distance):
            self.__read_data_file(distance)
            print('Data file existed, skip for distance =', distance)
            return
        assert(distance > 1)
        for permutation in self.distance_map[distance-1]:
            for reversal in self.__reversal_generator(self.__unhash(permutation)):
                if self.__contains(reversal, distance-2):
                    continue
                self.__insert(reversal, distance)
        self.__write_data_file(distance)
    
    def __reversal_generator(self, x):
        for i, j in itertools.combinations(range(self.n), 2):
            yield x[:i] + x[i:j+1][::-1] + x[j+1:]

    def __contains(self, permutation, distance):
        return self.__hash(permutation) in self.distance_map[distance]
            
    def __insert(self, permutation, distance):
        self.distance_map[distance].add(self.__hash(permutation))
    
    def __hash(self, permutation):
        return '|'.join(map(str, permutation))
    
    def __unhash(self, hashed):
        return hashed.split('|')
    
def main():
    query = ReversalDistanceQuery(10)
    
    #pairs = read_permutation_pairs()
    #normalized_permutations = [normalize(pair) for pair in pairs]
    #distances = []
    #for p in normalized_permutations:
    #    distances.append(breakpoint_reversal_distance(p))
    #print(' '.join(map(str, distances)))
    
if __name__ == '__main__':
    sys.exit(main())
