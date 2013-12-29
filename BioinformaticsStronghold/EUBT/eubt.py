import copy
import sys

class UnrootedBinaryTreeEnumerator():
    def __init__(self, species):
        self.species = species
        self.internal_node_id = 0
        self.enumerated_trees = []
        
    def enumerate(self):
        self.build_tree()
        for tree in self.enumerated_trees:
            yield self.newick_format(tree)
        
    def build_tree(self):
        tree = {}
        tree['Root'] = self.species[0:3]
        n = len(self.species)
        bfs_queue = [[tree, 3]]
        while bfs_queue:
            tree, num_leave = bfs_queue[0]
            bfs_queue.pop(0)
            if num_leave >= n:
                self.enumerated_trees.append(tree)
                continue
            for extended in self.extend_tree(tree, num_leave):
                next = [extended, num_leave + 1]
                bfs_queue.append(next)            
        
    def extend_tree(self, tree, num_leave):
        rv = []
        node = 'Root'
        dfs_stack = ['Root']
        while dfs_stack:
            u = dfs_stack[-1]
            dfs_stack.pop(-1)
            for v in tree[u]:
                next_tree = copy.deepcopy(tree)
                breaking_node = self.generate_internal_node()
                next_tree[u].remove(v)
                next_tree[u].append(breaking_node)
                next_tree[breaking_node] = [v, self.species[num_leave]]
                rv.append(next_tree)
                if v in tree:
                    dfs_stack.append(v)
        return rv
                
    def generate_internal_node(self):
        node = 'Internal[%d]' % self.internal_node_id
        self.internal_node_id += 1
        return node
    
    def newick_format(self, tree):
        return self.depth_first_search(tree, 'Root') + ';'
        
    def depth_first_search(self, tree, node):
        if node not in tree:
            return node
        rv = '('
        n = len(tree[node])
        for i in range(n):
            rv += self.depth_first_search(tree, tree[node][i])
            if i != n-1:
                rv += ','
        rv += ')'
        return rv
        
def main():
    species = sys.stdin.read().strip().split()
    enumerator = UnrootedBinaryTreeEnumerator(species)
    for tree in enumerator.enumerate():
        print(tree)
    
if __name__ == '__main__':
    sys.exit(main())
