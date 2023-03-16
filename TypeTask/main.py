from numpy import *
import scipy

# Выставляем количество знаков после запятой в матрице
set_printoptions(precision=3)

# Получаем вариант
N = int(input("Ваш вариант: "))

# Генерируем матрицу
MATRIX = matrix([
    [N, 4*N, N+2, N, N, N],
    [0, N+1, N, N, N+1, 2],
    [2*N, N, N+2, 2*N, 3*N, N+1],
    [N, 7*N, N+1, N, 4*N, N],
    [N, N+1, N, N, N+1, 2],
    [2*N, N, N+2, 2*N, N, N+1]
])

# Находим определитель матрицы
DET          = int(linalg.det(MATRIX))

MATRIX_UNION = f"{N} * {MATRIX / N}"

# Находим транспонированную матрицу
MATRIX_TRANSPOSE    = transpose(MATRIX)

# Находим сумму исходной матрицы и транспонированной
MATRIX_CONCAT       = MATRIX + MATRIX_TRANSPOSE

# Находим произведение исходной матрицы и транспонированной
MATRIX_MULTIPLY_D_T = MATRIX * MATRIX_TRANSPOSE

# Находим ранг матрицы
MATRIX_RANK         = linalg.matrix_rank(MATRIX)

# Находим инвертированную матрицу
MATRIX_INVERT       = linalg.inv(MATRIX)

# Находим произведение исходной матрицы и обратной
MATRIX_MULTIPLY_D_I = MATRIX * MATRIX_INVERT

# Находим собственные значения и вектора матрицы
MATRIX_EIG_W, MATRIX_EIG_V = linalg.eig(MATRIX)

# Проиизводим LU-разложение матрицы
P, L, U      = scipy.linalg.lu(MATRIX)
MATRIX_LU = L @ U

# print(MATRIX)
print(L)
print(U)
# print(MATRIX_INVERT)
print(MATRIX_LU)
