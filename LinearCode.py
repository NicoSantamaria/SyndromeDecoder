import MatrixMod


class LinearCode:
    def __init__(self, n: int = 2) -> None:
        self.mod = n
        self.G = None
        self.H = None
    

    def identity_matrix(self, n: int) -> [[int]]:
        """
        probably a neater way to do this
        """
        identity_mat = [[0 for i in range(n)] for j in range(n)]

        for i in range(n):
            identity_mat[i][i] = 1

        return identity_mat


    def generator_matrix(self, mat: [[int]]) -> None:
        GF = MatrixMod(mat, self.mod)
        GF.rref()
        self.G = GF.matrix

        num_rows = len(self.G)
        num_cols = len(self.G[0])
        dim = min(num_rows, num_cols)
        parity_dim = max(num_rows, num_cols) - dim

        non_identity = GF.get_transpose([GF.negate_vector(self.G[i][dim:]) for i in range(num_rows)])
        self.H = [non_identity[i] + self.identity_matrix(parity_dim)[i] for i in range(parity_dim)]


    def parity_check(self, mat: [[int]]) -> None:
        """
        Assumes that H is in reduced form
        """
        GF = MatrixMod(mat, self.mod)
        self.H = GF.matrix

        num_rows = len(self.H)
        num_cols = len(self.H[0])
        parity_dim = min(num_rows, num_cols)
        dim = max(num_rows, num_cols) - parity_dim

        non_identity = GF.get_transpose([GF.negate_vector(self.H[i][:dim]) for i in range(num_rows)])
        self.G = [self.identity_matrix(dim)[i] + non_identity[i] for i in range(dim)]


    def equivalent_codes(self, mat1: [[int]], mat2: [[int]]) -> bool:
        mat1_rref = MatrixMod(mat1, self.mod)
        mat2_rref = MatrixMod(mat2, self.mod)
        mat1_rref.rref()
        mat2_rref.rref()

        return mat1_rref.matrix == mat2_rref.matrix


    def get_syndrome(self) -> None:
        return
    

    def standard_array(self) -> None:
        return
    

# C = Code()
# C.generator_matrix([[1,0,0,0,1,0,1], [0,1,0,0,1,1,1], [0,0,1,0,1,1,0],[0,0,0,1,0,1,1]])
# C.parity_check([[1, 1, 1, 0, 1, 0, 0], [0, 1, 1, 1, 0, 1, 0], [1, 1, 0, 1, 0, 0, 1]])