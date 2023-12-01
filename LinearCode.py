"""
What remains:
    1. Handling exceptions
    2. Making code more readable/efficient/pretty
    3. Comments

Then it's ready to be implemented into the website
"""

from MatrixMod import *


class LinearCode:
    def __init__(self, matrix: [[int]]=None, n: int=2, is_generator: bool=True) -> None:
        self.mod = n
        self.GF = MatrixMod(self.mod, matrix)
        self.H = None
        self.G = None
        
        self.generator_matrix() if is_generator else self.parity_check()

        num_rows = len(self.G[0])
        self.coset_leaders = self.get_coset_leaders(num_rows, self.mod - 1)
        self.syndromes = [self.get_syndrome(vector) for vector in self.coset_leaders]
    

    def identity_matrix(self, n: int, multiple: int=1) -> [[int]]:
        return  [[0 if i != j else multiple for i in range(n)] for j in range(n)]
    

    def decode(self, vector: [int]) -> [int]:
        syndrome = self.get_syndrome(vector)
        pos = self.syndromes.index(syndrome)
        leader = self.coset_leaders[pos]

        return self.GF.reduce_vector([v - l for v, l in zip(vector, leader)])


    def generator_matrix(self) -> None:
        self.GF.rref()
        self.G = self.GF.matrix

        num_rows = len(self.G)
        num_cols = len(self.G[0])
        dim = min(num_rows, num_cols)
        parity_dim = max(num_rows, num_cols) - dim

        non_identity = self.GF.get_transpose([
            self.GF.negate_vector(self.G[i][dim:]) 
            for i in range(num_rows)])
        self.H = [non_identity[i] + self.identity_matrix(parity_dim)[i] for i in range(parity_dim)]


    def parity_check(self) -> None:
        """
        Assumes that H is in reduced form
        """
        self.H = self.GF.matrix

        num_rows = len(self.H)
        num_cols = len(self.H[0])
        parity_dim = min(num_rows, num_cols)
        dim = max(num_rows, num_cols) - parity_dim

        non_identity = self.GF.get_transpose([
            self.GF.negate_vector(self.H[i][:dim]) 
            for i in range(num_rows)])
        self.G = [self.identity_matrix(dim)[i] + non_identity[i] for i in range(dim)]


    def get_syndrome(self, vector: [int]) -> [int]:
        return self.GF.reduce_vector([self.GF.dot_product(vector, col) for col in self.H])


    def get_coset_leaders(self, length: int, n: int) -> [[int]]:
        if n == 0:
            return [[0] * length]

        else:
            return self.identity_matrix(length, n) + self.get_coset_leaders(length, n-1)

# Ham(3,3)
Hamming = LinearCode([[0,0,0,0,1,1,1,1,1,1,1,1,1],
                      [0,1,1,1,0,0,0,1,1,1,2,2,2],
                      [1,0,1,2,0,1,2,0,1,2,0,1,2]], 3, False)