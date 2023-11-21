class MatrixMod:
    def __init__(self, m: [[int]], p: int = 0) -> None:
        """
        Note: for now, q must be prime
        """
        self.matrix = m
        self.modulus = p

    def __str__(self) -> str:
        return

    def reduce_vector(self, vector: [int]) -> [int]:
        return [vector[i] % self.modulus for i in range(len(vector))]

    def multiply(self, x: int, y: int) -> int:
        return (x * y) % self.modulus
    
    def negate_vector(self, vector: [int]) -> [int]:
        return self.reduce_vector([self.multiply(-1, vector[i]) for i in range(len(vector))])

    def find_mult_inverse(self, num: int) -> int:
        """
        We assume the inverse of num exists since self.modulus is prime

        Can be optimized using Euler's Division Algorithm
        """
        for x in range(1, self.modulus):
            if self.multiply(x, num) == 1:
                return x

    def vector_addition(self, vector1: [int], vector2: [int]) -> [int]:
        if len(vector1) == len(vector2):
            return self.reduce_vector([vector1[i] + vector2[i] for i in range(len(vector1))])

        else:
            print("These vectors cannot be added.")
    
    def rref(self) -> None:
        """
        row reduce self.matrix in-place
        """
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
    
    def get_transpose(self, m: [[int]]) -> [[int]]:
        num_rows = len(m)
        num_cols = len(m[0])

        return [[m[j][i] for j in range(num_rows)] for i in range(num_cols)]

    