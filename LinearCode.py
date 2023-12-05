"""
Notes:
Check to see that it doesn't assume minimal weight = 1

Implement sympy for matrix multiplication
"""
from MatrixMod import *


class LinearCode:
    def __init__(self, matrix: [[int]]=None, n: int=2, is_generator: bool=True) -> None:
        """
        matrix: [[int]] - generator or parity check matrix for linear code
        n: int - modulus of finite field
        is_generator: bool - True is matrix is a generator matrix, false if matrix is a parity check
        return: None
        """
        # Call MatrixMod for linear algebra operations over finite field
        self.GF = MatrixMod(n, matrix); self.mod = n
        
        # Build self.G and self.H according to type of matrix input
        if is_generator:
            self.generator_matrix()
            num_rows = len(self.G[0])
        else:
            self.H = matrix
            num_rows = len(self.H) + len(self.H[0])

        # Build and store coset leaders and syndromes
        self.coset_leaders = self.get_coset_leaders(num_rows, self.mod - 1)
        self.syndromes = [self.get_syndrome(vector) for vector in self.coset_leaders]

    

    def decode(self, vector: [int]) -> [int]:
        """
        vector: [int] - message to be decoded
        return: [int] - error-corrected code
        """
        # Find associated coset leader for message
        # The associated coset leader has the same syndrome as the message 
        syndrome = self.get_syndrome(vector)
        pos = self.syndromes.index(syndrome)
        leader = self.coset_leaders[pos]

        try:
            # Decode with (vector) - (coset leader)
            return self.GF.reduce_vector([v - l for v, l in zip(vector, leader)])
        except:
            print("Decode error.")

    
    def get_min_weight(self) -> int:
        return 

    
    def identity_matrix(self, n: int, multiple: int=1) -> [[int]]:
        """
        n: int - length of matrix
        multiple: int - integer multiple of matrix
        return: [[int]] - nxn identity matrix multiplied by multiple
        """
        return  [[0 if i != j else multiple for i in range(n)] for j in range(n)]


    def generator_matrix(self) -> None:
        """
        Modifies self.G in-place into standard form, and 
        builds self.H from self.H
        return: None
        """
        # Call MatrixMod to put self.G in standard form
        self.GF.rref()
        self.G = self.GF.matrix

        num_rows = len(self.G)
        num_cols = len(self.G[0])
        dim = min(num_rows, num_cols)
        parity_dim = max(num_rows, num_cols) - dim

        try:
            # Constructs the portion of self.H which is not
            # equivalent to the identity matrix
            non_identity = self.GF.get_transpose([
                self.GF.negate_vector(self.G[i][dim:])
                for i in range(num_rows)
                ])
            
            # Joins the non-identity portion with the appropriately sized
            # identity matrix to build self.H
            self.H = [
                non_identity[i] + self.identity_matrix(parity_dim)[i] 
                for i in range(parity_dim)
                ]
        except:
            print("Generator matrix error.")


    def get_syndrome(self, vector: [int]) -> [int]:
        """
        vector: [int] - a vector in the codespace
        return: [int] - syndrome associated with the vector
        """
        try:
            # The syndrome is equivalent to the vector multiplied
            # by the transpose of the parity check matrix
            return self.GF.reduce_vector([
                self.GF.dot_product(vector, col)
                for col in self.H
                ])
        except:
            print("Syndrome error.")
            return None


    def get_coset_leaders(self, length: int, n: int) -> [[int]]:
        """
        length: int - length of the codewords
        n: int - modulus of the field, distinct from self.mod for recursion
        return: [[int]] - list containing the codewords of minimal weight in the code
        """
        if n == 0:
            return [[0] * length]
        
        # For a linear code, minimum possible weight is 1
        return self.identity_matrix(length, n) + self.get_coset_leaders(length, n-1)


# Ham(3,3)
Hamming = LinearCode([[0,0,0,0,1,1,1,1,1,1,1,1,1],[0,1,1,1,0,0,0,1,1,1,2,2,2],[1,0,1,2,0,1,2,0,1,2,0,1,2]], 3, False)
print(Hamming.decode([1,3,4,4,1,2,3,4,5,5,1,3,2]))