import itertools
import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    n = int(lines[0])
    edge_list_array = []
    is_begin = True
    num_edge = None
    edge_list = None
    for line in lines[1:]:
        if not line:
            continue
        edge = list(map(int, line.split()))
        if is_begin:
            is_begin = False
            edge_list = [edge]
            num_edge = edge[1]
            continue    
        edge_list.append(edge)
        if len(edge_list) - 1 == num_edge:
            is_begin = True
            edge_list_array.append(edge_list)
    assert(n == len(edge_list_array))
    return edge_list_array
    
class Graph():
    def __init__(self, edge_list):
        self.edge_list = edge_list
        self.graph_impl = None
        self.__build_graph()
    
    def has_square(self):
        n = self.edge_list[0][0]
        aux_matrix = [[False for i in range(n+1)] for j in range(n+1)]
        for i in range(1, n+1):
            for j, k in itertools.combinations(self.graph_impl[i], 2):
                if aux_matrix[j][k] or aux_matrix[k][j]:
                    return True
                else:
                    aux_matrix[j][k] = True
                    aux_matrix[k][j] = True
        return False
        
    def __build_graph(self):
        n = self.edge_list[0][0]
        self.graph_impl = {}
        for node in range(n):
            self.graph_impl[node + 1] = []
        for x, y in self.edge_list[1:]:
            self.graph_impl[x].append(y)
            self.graph_impl[y].append(x)
        
def main():
    edge_list_array = read_dataset()
    rv = []
    for edge_list in edge_list_array:
        graph = Graph(edge_list)
        if graph.has_square():
            rv.append(1)
        else:
            rv.append(-1)
    print(' '.join(map(str, rv)))
    
if __name__ == '__main__':
    sys.exit(main())
