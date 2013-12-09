import sys

class LongestRepeatedSubstringFinder():
    def __init__(self, dna_string, k, edges):
        self.dna_string = dna_string
        self.k = k
        self.edges = edges
        
        self.tree = {}
        self.root_nodes = set()
        self.leaf_nodes = set()
        
        self.repeated_node_map = {}
        
    def find(self):
        self.__build_tree()
        self.__build_repeated_node_map()
        return self.__build_longest_repeated_substring()

    def __build_tree(self):
        internal_nodes = set()
        for parent_node, child_node, location, length in self.edges:
            self.tree[child_node] = [parent_node, location, length]
            self.root_nodes.add(parent_node)
            internal_nodes.add(child_node)
            self.leaf_nodes.add(child_node)
            if parent_node in self.leaf_nodes:
                self.leaf_nodes.remove(parent_node)
        internal_nodes -= self.leaf_nodes
        self.root_nodes -= internal_nodes
        assert(len(self.root_nodes) == 1)        
    
    def __build_repeated_node_map(self):
        for edge in self.tree:
            self.repeated_node_map[edge] = 0
            self.repeated_node_map[self.tree[edge][0]] = 0
        for leaf in self.leaf_nodes:
            curr_parent_node = leaf
            while curr_parent_node not in self.root_nodes:
                self.repeated_node_map[curr_parent_node] += 1
                curr_parent_node = self.tree[curr_parent_node][0]
            self.repeated_node_map[curr_parent_node] += 1
    
    def __build_longest_repeated_substring(self):
        longest_so_far = ''
        for node in self.repeated_node_map:
            if self.repeated_node_map[node] == self.k:
                substring = self.__build_substring(node) 
                if len(substring) > len(longest_so_far):
                    longest_so_far = substring
        return longest_so_far
        
    def __build_substring(self, node):
        curr_parent_node = node
        substring = ''
        while curr_parent_node not in self.root_nodes:
            edge_info = self.tree[curr_parent_node]
            i, j = edge_info[1] - 1, edge_info[1] + edge_info[2] - 1
            substring = self.dna_string[i:j] + substring
            curr_parent_node = edge_info[0]
        return substring
    
def read_dataset():
    dataset = sys.stdin.readlines()
    dna_string = dataset[0].strip()
    k = int(dataset[1].strip())
    edges = []
    for i in range(2, len(dataset)):
        parent_node, child_node, location, length = dataset[i].strip().split()
        location = int(location)
        length = int(length)
        edges.append([parent_node, child_node, location, length])
    return dna_string, k, edges    
            
def main():
    dna_string, k, edges = read_dataset()
    finder = LongestRepeatedSubstringFinder(dna_string, k, edges)
    print(finder.find())
    
if __name__ == '__main__':
    sys.exit(main())
