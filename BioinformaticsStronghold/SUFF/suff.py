import sys

class SuffixTree():
    def __init__(self, string):
        self.suffix_tree = {}
        self.__build(string)
        self.edges = []
        
    def get_edges(self):
        self.__build_edges_recursively(self.suffix_tree)
        return self.edges
        
    def __build(self, string):
        n = len(string)
        for i in range(n):
            self.__insert(string[i:])
        
    def __insert(self, suffix):
        curr_node = self.suffix_tree
        curr_suffix = suffix
        while True:
            matched_node, matched_lcp = self.__search_child_node(curr_suffix, curr_node)
            if not matched_lcp:
                curr_node[curr_suffix] = {}
                break
            n = len(matched_lcp)
            if n != len(matched_node):
                matched_node_children = curr_node[matched_node]            
                curr_node[matched_lcp] = {
                    curr_suffix[n:] : {}, 
                    matched_node[n:] : matched_node_children,
                }
                curr_node.pop(matched_node)
                break
            curr_node = curr_node[matched_node]
            curr_suffix = curr_suffix[n:]
            
    def __search_child_node(self, string, parent_node):
        for child_node in parent_node:
            lcp = self.__get_longest_common_prefix(string, child_node)
            if lcp:
                return child_node, lcp
        return None, None
            
    def __get_longest_common_prefix(self, x, y):
        n = min(len(x), len(y))
        longest_common_pos = n
        for i in range(n):
            if x[i] != y[i]:
                longest_common_pos = i
                break
        return x[:longest_common_pos]
    
    def __build_edges_recursively(self, parent_node): 
        for child_node in parent_node:
            self.edges.append(child_node)
            self.__build_edges_recursively(parent_node[child_node])
    
def main():
    dna_string = sys.stdin.readline().strip()
    tree = SuffixTree(dna_string)
    print('\n'.join(tree.get_edges()))
    
if __name__ == '__main__':
    sys.exit(main())
