class MatrixMod:
    def __init__(self, p: int, m: [[int]]= None) -> None:
        """
        Note: for now, q must be prime
        """
        self.matrix = m
        self.modulus = p


    def reduce_vector(self, vector: [int]) -> [int]:
        return [x % self.modulus for x in vector]


    def multiply(self, x: int, y: int) -> int:
        return (x * y) % self.modulus
    

    def negate_vector(self, vector: [int]) -> [int]:
        return self.reduce_vector([self.multiply(-1, vector[i]) for i in range(len(vector))])


    def dot_product(self, vector1: [int], vector2: [int]) -> int:
        return sum(self.multiply(a, b) for a, b in zip(vector1, vector2))


    def multiply_matrix(self, mat1: [[int]], mat2: [[int]]) -> [[int]]:
        return [[self.dot_product(row, col) for row in mat1] for col in self.get_transpose(mat2)]

    # Can be optimized for large num with Euler's Division Algorithm
    def find_mult_inverse(self, num: int) -> int:
        for x in range(1, self.modulus):
            if self.multiply(x, num) == 1:
                return x


    def vector_addition(self, vector1: [int], vector2: [int]) -> [int]:
        return self.reduce_vector([v1 + v2 for v1, v2 in zip(vector1, vector2)])
    
    # Row reduce self.matrix in place
    def rref(self) -> None:
        num_rows = len(self.matrix)
        num_cols = len(self.matrix[0])

        i, j = 0, 0,

        while True:
            if i >= num_rows or j >= num_cols:
                break
                
            if self.matrix[i][j] == 0:
                non_zero_row = i

                while non_zero_row < num_rows and self.matrix[non_zero_row][j] == 0:
                    non_zero_row += 1
                
                if non_zero_row == num_rows:
                    j += 1
                    continue
            
                self.matrix[i], self.matrix[non_zero_row] = self.matrix[non_zero_row], self.matrix[i]

            pivot = self.matrix[i][j]
            pivot_inverse = self.find_mult_inverse(pivot)
            self.matrix[i] = [self.multiply(x, pivot_inverse) for x in self.matrix[i]]

            for other_row in range(num_rows):
                if other_row == i:
                    continue

                if self.matrix[other_row][j] != 0:
                    self.matrix[other_row] = self.reduce_vector([y - self.matrix[other_row][j]*x
                                              for (x,y) in zip(self.matrix[i], self.matrix[other_row])])
            
            i += 1; j += 1
    

    def get_transpose(self, mat: [[int]]) -> [[int]]:
        num_rows = len(mat)
        num_cols = len(mat[0])

        return [[mat[j][i] for j in range(num_rows)] for i in range(num_cols)]
    
# M = MatrixMod(5)
