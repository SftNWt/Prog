from numpy import *
import scipy

N = int(input("Ваш вариант: "))

MATRIX = array([
    [N, 4*N, N+2, N, N, N],
    [0, N+1, N, N, N+1, 2],
    [2*N, N, N+2, 2*N, 3*N, N+1],
    [N, 7*N, N+1, N, 4*N, N],
    [N, N+1, N, N, N+1, 2],
    [2*N, N, N+2, 2*N, N, N+1]
])

set_printoptions(precision=3)


DET          = int(linalg.det(MATRIX))
MATRIX_UNION = f"{N} * {MATRIX / N}"
MATRIX_TRANSPOSE    = transpose(MATRIX)
MATRIX_CONCAT       = MATRIX + MATRIX_TRANSPOSE
MATRIX_MULTIPLY_D_T = MATRIX * MATRIX_TRANSPOSE
MATRIX_RANK         = linalg.matrix_rank(MATRIX)
MATRIX_INVERT       = linalg.inv(MATRIX)
MATRIX_MULTIPLY_D_I = MATRIX * MATRIX_INVERT
MATRIX_EIG_W, MATRIX_EIG_V = linalg.eig(MATRIX)
P, MATRIX_L, MATRIX_U      = scipy.linalg.lu(MATRIX)
# MATRIX_LU           = concatenate((MATRIX_L, MATRIX_U))

print(MATRIX)
print(MATRIX_L * MATRIX_U * P)
