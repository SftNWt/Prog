from numpy import *
from functions import *

# Выставляем количество знаков после запятой в матрице
set_printoptions(precision=3)

# Получаем вариант
N = int(input("\n\033[34m\033[1mВаш вариант: "))

MATRIX_A = array([
    [10 * N, N, N + 2],
    [N, 9 * N, 3 * N],
    [2 * N, 3 * N, 11 * N],
])
MATRIX_B = array([
    [2 * N], [N], [3 * N]
])

Gauss_Method(MATRIX_A, MATRIX_B, output = True)
Crout_Method(MATRIX_A, MATRIX_B, output = True)
jacobi(MATRIX_A, MATRIX_B, iteration = 2, output = True)
# GAUSS_SAIDEL  = gauss_saidel(MATRIX_A, MATRIX_B, iteration=2, accuracy=1e-6)
