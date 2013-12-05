import itertools
import sys

class RestrictionMap():
    def __init__(self, L):
        self.L = L
        self.X = None
        
    def create(self):
        L = self.L
        L.sort()
        width = L[-1]
        X = [0, width]
        L = self.L[:-1]
        self.place(L, X, width)
    
    def place(self, L, X, width):
        # Found one solution is good enough
        if self.X:
            return
            
        # Found
        if not L:
            X.sort()
            self.X = X
            return

        L.sort()
        y = L[-1]
        left_delta = self.__minkowski_delta(y, X)
        if self.__multiset_contained(left_delta, L):
            X.append(y)
            L = self.__multiset_minus(L, left_delta)
            self.place(list(L), list(X), width)
            L = L + left_delta
            X.remove(y)
        right_delta = self.__minkowski_delta(width - y, X)
        if self.__multiset_contained(right_delta, L):
            X.append(width - y)
            L = self.__multiset_minus(L, right_delta)
            self.place(list(L), list(X), width)
            L = L + right_delta
            X.remove(width - y)
        
    def __minkowski_delta(self, y, X):
        return [abs(y-x) for x in X]        
        
    def __multiset_minus(self, X, Y):
        for y in Y:
            X.remove(y)
        return X
    
    def __multiset_contained(self, X, Y):
        Z = list(Y)
        for x in X:
            if x not in Z:
                return False
            Z.remove(x)
        return True    

def main():
    difference_multiset = list(map(int, sys.stdin.read().strip().split()))
    m = RestrictionMap(difference_multiset)
    m.create()
    print(' '.join(map(str, m.X)))
        
if __name__ == '__main__':
    sys.exit(main())
