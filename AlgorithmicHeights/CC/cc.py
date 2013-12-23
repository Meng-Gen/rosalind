import sys

def read_dataset():
    return [list(map(int, line.strip().split())) for line in sys.stdin.readlines()]
    
class Graph():
    def __init__(self, edge_list):
        self.edge_list = edge_list
        self.num_nodes = edge_list[0][0]
        self.graph_impl = {}
        self.__build_graph()
    
    def num_connected_components(self):
        count = 0
        visited_nodes = set()
        for node in range(1, self.num_nodes + 1):
            if node in visited_nodes:
                continue
            dfs_stack = [node]
            while dfs_stack:
                top = dfs_stack[-1]
                dfs_stack.pop(-1)
                if top in visited_nodes:
                    continue
                visited_nodes.add(top)
                for next in self.graph_impl[top]:
                    dfs_stack.append(next)
            count += 1
        return count
        
    def __build_graph(self):
        for node in range(self.num_nodes):
            self.graph_impl[node + 1] = []
        for x, y in self.edge_list[1:]:
            self.graph_impl[x].append(y)
            self.graph_impl[y].append(x)
        
def main():
    edge_list = read_dataset()
    graph = Graph(edge_list)
    print(graph.num_connected_components())
    
if __name__ == '__main__':
    sys.exit(main())
