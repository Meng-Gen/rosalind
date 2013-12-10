import sys

class SuffixTree():
    def __init__(self, string):
        self.suffix_tree = ['', []]
        self.edges = []
        self.__build(string)
        
    def get_edges(self):
        self.__get_edges_recursively(self.suffix_tree, '')
        return [_ for _ in self.edges if _ != '']
        
    def __build(self, string):
        for i in range(len(string)):
            self.__insert(string[i:])

        
    def __insert(self, suffix):
        curr_node = self.suffix_tree
        for child in suffix:
            found_child = False
            for next_node in curr_node[1]:
                if next_node[0] == child:
                    found_child = True
                    break
            if not found_child:
                curr_node[1].append([child, []])
            for next_node in curr_node[1]:
                if next_node[0] == child:
                    curr_node = next_node
                    break
    
    def __get_edges_recursively(self, curr_node, queued_substring):
        queued_substring += curr_node[0]
        if len(curr_node[1]) == 0:
            return self.edges.append(queued_substring)
        elif len(curr_node[1]) == 1:
            return self.__get_edges_recursively(curr_node[1][0], queued_substring)
        self.edges.append(queued_substring)
        queued_substring = ''
        for next_node in curr_node[1]:
            self.__get_edges_recursively(next_node, queued_substring)
    
def main():
    dna_string = sys.stdin.readline().strip()
    tree = SuffixTree(dna_string)
    print('\n'.join(tree.get_edges()))
    
if __name__ == '__main__':
    sys.exit(main())
