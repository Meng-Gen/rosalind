import sys

class SuffixTree():
    def __init__(self, string):
        self.TERMINAL_SYMBOL = '$'
        self.string = string + self.TERMINAL_SYMBOL
        self.suffix_tree = {}
        self.maximal_repeats = None
        self.build_suffix_tree()
        
    def get_maximal_repeats(self):
        self.maximal_repeats = set()
        self.depth_first_search()
        return self.maximal_repeats
        
    def build_suffix_tree(self):
        n = len(self.string)
        self.insert_suffix(self.string, None)
        for i in range(1, n):
            self.insert_suffix(self.string[i:], self.string[i-1])
        
    def insert_suffix(self, suffix, left_character):
        curr_node = self.suffix_tree
        curr_suffix = suffix
        while True:
            matched_node, matched_lcp = self.search_child_node(curr_suffix, curr_node)
            if not matched_lcp:
                curr_node[curr_suffix] = [{}, left_character]
                break
            n = len(matched_lcp)
            if n != len(matched_node):
                matched_node_children = curr_node[matched_node]
                curr_node[matched_lcp] = {
                    curr_suffix[n:] : [{}, left_character], 
                    matched_node[n:] : matched_node_children,
                }
                curr_node.pop(matched_node)
                break
            curr_node = curr_node[matched_node]
            curr_suffix = curr_suffix[n:]
            
    def search_child_node(self, string, parent_node):
        for child_node in parent_node:
            lcp = self.get_longest_common_prefix(string, child_node)
            if lcp:
                return child_node, lcp
        return None, None
            
    def get_longest_common_prefix(self, x, y):
        n = min(len(x), len(y))
        longest_common_pos = n
        for i in range(n):
            if x[i] != y[i]:
                longest_common_pos = i
                break
        return x[:longest_common_pos]
    
    def depth_first_search(self): 
        root = [None, self.suffix_tree, []]
        dfs_stack = [root]
        while dfs_stack:
            u, children, path = dfs_stack[-1]
            dfs_stack.pop(-1)
            for v in children:
                next_path = path + [u]
                if isinstance(children, list):
                    continue
                next = [v, children[v], next_path]
                dfs_stack.append(next)
            if self.is_left_diverse(u, children):
                self.insert_maximal_repeat(''.join(path[1:] + [u]))
            
    def is_left_diverse(self, node, children):
        if self.has_internal_child(children):
            return False
        if node[-1] == self.TERMINAL_SYMBOL:
            return False
        left_characters = set()
        for v in children:
            left_characters.add(children[v][1])
        if len(left_characters) == 1:
            return False
        else:
            return True
   
    def has_internal_child(self, children):
        if isinstance(children, list):
            return False
        for v in children:
            if v[-1] != self.TERMINAL_SYMBOL:
                return True
        return False
    
    def insert_maximal_repeat(self, maximal_repeat):
        if len(maximal_repeat) >= 20:
            self.maximal_repeats.add(maximal_repeat)
   
def main():
    dna_string = sys.stdin.readline().strip()
    rv = set()
    string_stack = [dna_string]
    while string_stack:
        top = string_stack[-1]
        string_stack.pop(-1)
        tree = SuffixTree(top)
        maximal_repeats = tree.get_maximal_repeats()
        rv |= maximal_repeats
        string_stack += list(maximal_repeats)
    print('\n'.join(list(rv)))
    
if __name__ == '__main__':
    sys.exit(main())
