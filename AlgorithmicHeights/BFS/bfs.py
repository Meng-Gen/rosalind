import sys

def read_dataset():
    return [list(map(int, line.strip().split())) for line in sys.stdin.readlines()]
    
class Graph():
    def __init__(self, edge_list):
        self.edge_list = edge_list
        self.num_nodes = edge_list[0][0]
        self.graph_impl = {}
        self.__build_graph()
    
    def get_degree_array(self):
        return [len(self.graph_impl[x + 1]) for x in range(self.num_nodes)]

    def breadth_first_search(self):
        distance_array = [-1 for _ in range(self.num_nodes)]
        init_state = [1, 0]
        visited_nodes = set()
        bfs_queue = [init_state]
        while bfs_queue:
            top_node, top_distance = bfs_queue[0]
            bfs_queue = bfs_queue[1:]
            if top_node in visited_nodes:
                continue
            visited_nodes.add(top_node)
            distance_array[top_node - 1] = top_distance # 0-based
            for next_node in self.graph_impl[top_node]:
                next_distance = top_distance + 1
                bfs_queue.append([next_node, next_distance])
        return distance_array
        
    def __build_graph(self):
        for node in range(self.num_nodes):
            self.graph_impl[node + 1] = []
        for x, y in self.edge_list[1:]:
            self.graph_impl[x].append(y)
        
def main():
    edge_list = read_dataset()
    graph = Graph(edge_list)
    print(' '.join(map(str, graph.breadth_first_search())))
    
if __name__ == '__main__':
    sys.exit(main())
