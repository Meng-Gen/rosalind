import sys
    
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
            elif self.c0 and (self.c0.isalpha() or self.c0 == '_'):
                self.__scan_string()
            elif not self.c0:
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
        # Branch :: Subtree
        return ['Branch', self.__parse_subtree()]
        
    def __parse_name(self):
        # Name :: empty | string
        name = None
        token = self.scanner.get_current_token()
        if token[0] == 'TK_STRING':
            name = token[1]
            self.__consume_token()
        return ['Name', name]
    
    def __consume_token(self):
        self.scanner.move_to_next_token()
        
class Pedigree():
    def __init__(self, newick_format):
        self.abstract_syntax_tree = None
        self.curr_id = 0
        self.node_name_map = {}
        self.edges = {}
        self.__build(newick_format)
    
    def get_genotype_prob(self):
        return self.__get_genotype_prob_recursively(0)
    
    def __build(self, newick_format):
        parser = NewickFormatParser(newick_format)
        self.abstract_syntax_tree = parser.parse()
        self.__interpret(self.abstract_syntax_tree)
        
    def __interpret(self, block): 
        self.__interpret_subtree(self.abstract_syntax_tree[1], None)
        
    def __interpret_subtree(self, subtree_block, parent_id):
        type = subtree_block[1][0]
        if type == 'Internal':
            self.__interpret_internal(subtree_block[1], parent_id)
        elif type == 'Leaf':
            self.__interpret_leaf(subtree_block[1], parent_id)
            
    def __interpret_internal(self, internal_block, parent_id): 
        name, next_id = internal_block[2][1], self.curr_id
        self.__update_graph(name, parent_id)
        self.__interpret_branch_list(internal_block[1], next_id)
        
    def __interpret_branch_list(self, branch_list_block, parent_id): 
        self.__interpret_branch(branch_list_block[1], parent_id)
        if len(branch_list_block) == 3:
            self.__interpret_branch_list(branch_list_block[2], parent_id)
        
    def __interpret_branch(self, branch_block, parent_id): 
        self.__interpret_subtree(branch_block[1], parent_id)
        
    def __interpret_leaf(self, leaf_block, parent_id): 
        name = leaf_block[1][1]
        self.__update_graph(name, parent_id)
        
    def __update_graph(self, name, parent_id):
        self.node_name_map[self.curr_id] = name
        if parent_id is not None:
            if parent_id not in self.edges:
                self.edges[parent_id] = set()
            if self.curr_id not in self.edges:
                self.edges[self.curr_id] = set()
            self.edges[parent_id].add(self.curr_id)
        self.curr_id += 1
    
    def __get_genotype_prob_recursively(self, parent_id):
        children = self.edges[parent_id]
        if not children:
            genotype = self.node_name_map[parent_id]
            if genotype == 'AA':
                return [1.0, 0.0, 0.0]
            elif genotype == 'Aa':
                return [0.0, 1.0, 0.0]
            elif genotype == 'aa':
                return [0.0, 0.0, 1.0]
        left, right = children
        x = self.__get_genotype_prob_recursively(left)
        y = self.__get_genotype_prob_recursively(right)
        total_prob = [
            x[0]*y[0] + x[0]*y[1]/2 + x[1]*y[0]/2 + x[1]*y[1]/4,
            x[0]*y[1]/2 + x[1]*y[0]/2 + x[2]*y[1]/2 + x[1]*y[2]/2 + x[2]*y[0] + x[0]*y[2] + x[1]*y[1]/2,
            x[2]*y[2] + x[2]*y[1]/2 + x[1]*y[2]/2 + x[1]*y[1]/4,
        ]
        return total_prob
    
def main():
    newick_format = sys.stdin.readline().strip()
    pedigree = Pedigree(newick_format)
    print(' '.join(map(str, pedigree.get_genotype_prob())))
    
if __name__ == '__main__':
    sys.exit(main())
