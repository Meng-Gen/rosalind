import sys

def read_dataset():
    lines = [_.strip() for _ in sys.stdin.readlines()]
    problem_array = []
    curr_pos = 0
    while curr_pos < len(lines):
        if lines[curr_pos]:
            newick_format = lines[curr_pos]
            edge = lines[curr_pos+1].split()
            problem_array.append([newick_format, edge])
            curr_pos += 2
        else:
            curr_pos += 1
    return problem_array
    
class NewickFormatScanner():
    def __init__(self, source):
        self.source = source
        self.pos = -1
        self.c0 = None
        self.tokens = []
        self.curr_token_pos = 0
    
    def scan(self):
        while True:
            self.__advance()
            if self.c0 == ';':
                self.tokens.append(['TK_SEMICOLON'])
            elif self.c0 == '(':
                self.tokens.append(['TK_LPAREN'])
            elif self.c0 == ')':
                self.tokens.append(['TK_RPAREN'])
            elif self.c0 == ',':
                self.tokens.append(['TK_COMMA'])
            elif self.c0 == ':':
                self.tokens.append(['TK_COLON'])
            elif self.c0:
                if self.c0.isalpha() or self.c0 == '_':
                    self.__scan_string()
                elif self.c0.isdigit():
                    self.__scan_number()
            if not self.c0:
                break
    
    def get_current_token(self):
        return self.tokens[self.curr_token_pos]

    def move_to_next_token(self):
        self.curr_token_pos += 1
        
    def __scan_string(self):
        string = self.c0
        while True:
            self.__advance()
            if self.c0 and (self.c0.isalpha() or self.c0 == '_'):
                string += self.c0
            else:
                self.__pushback()
                break
        self.tokens.append(['TK_STRING', string])
        
    def __scan_number(self):
        number = self.c0
        while True:
            self.__advance()
            if self.c0 and self.c0.isdigit():
                number += self.c0
            else:
                self.__pushback()
                break
        self.tokens.append(['TK_NUMBER', int(number)])    
        
    def __advance(self):
        self.pos += 1
        if self.pos < len(self.source):
            self.c0 = self.source[self.pos]
        else:
            self.c0 = None
            
    def __pushback(self):
        self.pos -= 1
        if self.pos >= 0:
            self.c0 = self.source[self.pos]
        else:
            self.c0 = None
            
class NewickFormatParser():
    def __init__(self, source):
        self.scanner = NewickFormatScanner(source)
        
    def parse(self):
        self.scanner.scan()
        return self.__parse_tree()
         
    def __parse_tree(self): 
        # Tree :: Subtree ";"
        subtree = self.__parse_subtree()
        self.__consume_token()
        return ['Tree', subtree]
        
    def __parse_subtree(self):
        # Subtree :: Leaf | Internal
        if self.scanner.get_current_token()[0] == 'TK_LPAREN':
            return ['Subtree', self.__parse_internal()]
        else:
            return ['Subtree', self.__parse_leaf()]
        
    def __parse_leaf(self):
        # Leaf :: Name
        return ['Leaf', self.__parse_name()]
        
    def __parse_internal(self):
        # Internal :: "(" BranchList ")" Name
        self.__consume_token()
        branch_list = self.__parse_branch_list()
        self.__consume_token()
        name = self.__parse_name()
        return ['Internal', branch_list, name]
        
    def __parse_branch_list(self):
        # BranchList :: Branch | Branch "," BranchList
        branch = self.__parse_branch()
        token = self.scanner.get_current_token()
        if token[0] == 'TK_COMMA':
            self.__consume_token()
            branch_list = self.__parse_branch_list()
            return ['BranchList', branch, branch_list]
        else:
            return ['BranchList', branch]
        
    def __parse_branch(self):
        # Branch :: Subtree Length
        subtree = self.__parse_subtree()
        length = self.__parse_length()
        return ['Branch', subtree, length]
        
    def __parse_name(self):
        # Name :: empty | string
        name = None
        token = self.scanner.get_current_token()
        if token[0] == 'TK_STRING':
            name = token[1]
            self.__consume_token()
        return ['Name', name]
    
    def __parse_length(self):
        # Length :: empty | ":" number
        length = None
        token = self.scanner.get_current_token()
        if token[0] == 'TK_COLON':
            self.__consume_token()
            number_token = self.scanner.get_current_token()
            length = number_token[1]
            self.__consume_token()
        return ['Length', length]
    
    def __consume_token(self):
        self.scanner.move_to_next_token()
        
class NewickFormatTree():
    def __init__(self, newick_format):
        self.abstract_syntax_tree = None
        self.curr_id = 0
        self.node_name_map = {}
        self.node_name_inversed_map = {}
        self.edges = {}
        self.edge_weight_map = {}
        self.__build(newick_format)
    
    def get_distance(self, edge):
        u, v = [self.node_name_inversed_map[_] for _ in edge]
        init_state = [u, 0]
        visited_nodes = set()
        bfs_queue = [init_state]
        while bfs_queue:
            top_id, top_distance = bfs_queue[0]
            bfs_queue = bfs_queue[1:]
            if top_id == v:
                return top_distance
            if top_id in visited_nodes:
                continue
            visited_nodes.add(top_id)
            for next_id in self.edges[top_id]:
                next_distance = top_distance + self.__get_weight(top_id, next_id)
                bfs_queue.append([next_id, next_distance])
        return None
    
    def __build(self, newick_format):
        parser = NewickFormatParser(newick_format)
        self.abstract_syntax_tree = parser.parse()
        self.__interpret(self.abstract_syntax_tree)
        
    def __interpret(self, block): 
        self.__interpret_subtree(self.abstract_syntax_tree[1], None, None)
        
    def __interpret_subtree(self, subtree_block, parent_id, length):
        type = subtree_block[1][0]
        if type == 'Internal':
            self.__interpret_internal(subtree_block[1], parent_id, length)
        elif type == 'Leaf':
            self.__interpret_leaf(subtree_block[1], parent_id, length)
            
    def __interpret_internal(self, internal_block, parent_id, length): 
        name, next_id = internal_block[2][1], self.curr_id
        self.__update_graph(name, parent_id, length)
        self.__interpret_branch_list(internal_block[1], next_id)
        
    def __interpret_branch_list(self, branch_list_block, parent_id): 
        self.__interpret_branch(branch_list_block[1], parent_id)
        if len(branch_list_block) == 3:
            self.__interpret_branch_list(branch_list_block[2], parent_id)
        
    def __interpret_branch(self, branch_block, parent_id): 
        length = branch_block[2][1]
        self.__interpret_subtree(branch_block[1], parent_id, length)
        
    def __interpret_leaf(self, leaf_block, parent_id, length): 
        name = leaf_block[1][1]
        self.__update_graph(name, parent_id, length)
        
    def __update_graph(self, name, parent_id, length):
        u, v = parent_id, self.curr_id
        self.node_name_map[v] = name
        if name:
            self.node_name_inversed_map[name] = v
        if u is not None:
            if u not in self.edges:
                self.edges[u] = set()
            if v not in self.edges:
                self.edges[v] = set()
            self.edges[u].add(v)
            self.edges[v].add(u)
            self.edge_weight_map[self.__get_edge_hash(u, v)] = length
            self.edge_weight_map[self.__get_edge_hash(v, u)] = length
        self.curr_id += 1
    
    def __get_edge_hash(self, u, v):
        return '''%d-%d''' % (u, v)
        
    def __get_weight(self, u, v):
        return self.edge_weight_map[self.__get_edge_hash(u, v)]
        
def main():
    distance_array = []
    for problem in read_dataset():
        newick_format, edge = problem
        tree = NewickFormatTree(newick_format)
        distance_array.append(tree.get_distance(edge))
    print(' '.join(map(str, distance_array)))
    
if __name__ == '__main__':
    sys.exit(main())
