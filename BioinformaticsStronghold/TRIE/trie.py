import sys

class Trie():
    def __init__(self):
        self.curr_id = 1
        self.impl = ['root', []]
        
    def insert(self, s):
        curr_node = self.impl
        for child in s:
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

    def print_adjacency_list(self):
        self.curr_id = 1
        self._print_adjacency_list_recursively(self.impl, 1)

    def _print_adjacency_list_recursively(self, curr_node, parent_id):
        for next_node in curr_node[1]:
            print(parent_id, self.curr_id + 1, next_node[0])
            self.curr_id += 1
            self._print_adjacency_list_recursively(next_node, self.curr_id)                    
    
def main():
    t = Trie()
    for dna in sys.stdin.readlines():
        t.insert(dna.strip())
    t.print_adjacency_list()
    
if __name__ == '__main__':
    sys.exit(main())
