from numpy import *
import math
import scipy

def forward_elimination(Matrix_A, Matrix_B, matrix_shape, output):

    # ?------------------------ Calculation Block ------------------------
    for row in range(0, matrix_shape - 1):
        for i in range(row + 1, matrix_shape):
            factor = Matrix_A[i, row] / Matrix_A[row, row]
            for j in range(row, matrix_shape):
                Matrix_A[i, j] = Matrix_A[i, j] - factor * Matrix_A[row, j]

            Matrix_B[i] = Matrix_B[i] - factor * Matrix_B[row]
    # ?------------------------ Calculation Block ------------------------

    # ?-------------------------- Output Block ---------------------------
        if output:
            string_matrix = 'Equation = '
            string = ''
            for i in range(0, len(Matrix_A)):
                if i == 0 or i == 2:
                    string += '\033[0m           '
                    string += f'\033[33m\033[1m[ {str(Matrix_A[i])[1:-1]}'
                    string += f' | x{i + 1} ]\033[0m   '
                    string += f'\033[33m\033[1m{str(Matrix_B[i])}'
                    string += '\n'
                elif i == 1:
                    string += string_matrix
                    string += f'\033[33m\033[1m[ {str(Matrix_A[i])[1:-1]}'
                    string += f' | x{i + 1} ]\033[33m\033[1m = '
                    string += f'\033[33m\033[1m{str(Matrix_B[i])}'
                    string += '\n'
            print(string)
    # ?-------------------------- Output Block ---------------------------

    return Matrix_A, Matrix_B


def back_substitution(Matrix_A, Matrix_B, matrix_shape):
    # ?------------------------ Calculation Block ------------------------
    matrix_zeros = zeros((matrix_shape, 1))
    matrix_zeros[matrix_shape - 1] = Matrix_B[matrix_shape - 1] / Matrix_A[matrix_shape - 1, matrix_shape - 1]
    for row in range(matrix_shape - 2, -1, -1):
        sums = Matrix_B[row]
        for j in range(row + 1, matrix_shape):
            sums = sums - Matrix_A[row, j] * matrix_zeros[j]
        matrix_zeros[row] = sums / Matrix_A[row, row]
    # ?------------------------ Calculation Block ------------------------

    return matrix_zeros


def Gauss_Method(Matrix_A, Matrix_B, output = False):
    if output:
        print('\n\033[31m\033[1m\033[3m------------ Start Gauss method -------------')
        print('\n')


    # ?------------------------ Calculation Block ------------------------
    matrix_shape = Matrix_A.shape[0]
    if any(diag(Matrix_A) == 0):
        raise ZeroDivisionError(('Division by zero will occur; '
                                 'pivoting currently not supported'))

    Matrix_A, Matrix_B = forward_elimination(Matrix_A, Matrix_B, matrix_shape, output)
    roots = back_substitution(Matrix_A, Matrix_B, matrix_shape)
    # ?------------------------ Calculation Block ------------------------



    if output:
    # ?-------------------------- Output Block ---------------------------
        string_matrix = 'Answer   = '
        string = ''
        rootsVar = ['x1', 'x2', 'x3']
        for i in range(0, len(Matrix_A)):
            if i == 0 or i == 2:
                string += '\033[0m           '
                string += f'\033[32m\033[1m[ {rootsVar[i]} = {"%.5f" % roots[i]} ]\n'
            elif i == 1:
                string += string_matrix
                string += f'\033[32m\033[1m[ {rootsVar[i]} = {"%.5f" % roots[i]} ]\n'
        print(string)
        print('\n\033[31m\033[1m\033[3m------------- End Gauss method --------------')
    # ?-------------------------- Output Block ---------------------------
    else:
        return roots


def Crout_Method(Matrix_A, Matrix_B, output = False):
    if output:
        print('\n\033[31m\033[1m\033[3m------------ Start Crout method -------------')
        print('\n')

    # ?------------------------ Calculation Block ------------------------
    cout = 0
    Count_Rows_Matrix, Count_Columns_Matrix = Matrix_A.shape
    if (Count_Rows_Matrix != Count_Columns_Matrix):
        if output:
            print('The matrix is not homogeneous')
        else:
            return 0

    L_Matrix = zeros((Count_Columns_Matrix, Count_Columns_Matrix))
    U_Matrix = zeros((Count_Columns_Matrix, Count_Columns_Matrix))
    Time_Matrix = 0

    for i in range(Count_Columns_Matrix):
        L_Matrix[i][0] = Matrix_A[i][0]
        U_Matrix[i][i] = 1
    for j in range(1, Count_Columns_Matrix):
        U_Matrix[0][j] = Matrix_A[0][j] / L_Matrix[0][0]
    for k in range(0, Count_Columns_Matrix):
        for i in range(k, Count_Columns_Matrix):
            for r in range(k):
                Time_Matrix += L_Matrix[i][r] * U_Matrix[r][k]
            L_Matrix[i][k] = Matrix_A[i][k] - Time_Matrix
            Time_Matrix = 0
        for j in range(k + 1, Count_Columns_Matrix):
            for r in range(k):
                Time_Matrix += L_Matrix[k][r] * U_Matrix[r][j]
            U_Matrix[k][j] = (Matrix_A[k][j] - Time_Matrix) / L_Matrix[k][k]
            Time_Matrix = 0


    y = zeros(Count_Columns_Matrix)
    y[0] = Matrix_B[0] / L_Matrix[0][0]
    for k in range(1, Count_Columns_Matrix):
        for r in range(k):
            Time_Matrix += L_Matrix[k][r] * y[r]
        y[k] = (Matrix_B[k] - Time_Matrix) / L_Matrix[k][k]
        print(y)
        Time_Matrix = 0

    x = zeros(Count_Columns_Matrix)
    x[Count_Columns_Matrix - 1] = y[Count_Columns_Matrix - 1]
    print(x)
    for k in range(Count_Columns_Matrix - 2, -1, -1):
        for r in range(k + 1, Count_Columns_Matrix):
            Time_Matrix += U_Matrix[k][r] * x[r]
        x[k] = y[k] - Time_Matrix
        Time_Matrix = 0
    print(x)

    roots = {}
    for i in range(Count_Columns_Matrix):
        roots["x" + str(i + 1)] = float("%.5f" % x[i])
    # print(roots)
    # ?------------------------ Calculation Block ------------------------


    print('\n\033[31m\033[1m\033[3m------------- End Crout method --------------')
    return [roots, L_Matrix, U_Matrix]


def jacobi(Matrix_A, Matrix_B, iteration = 2, Zeros_Matrix = None, output = False):
    if output:
        print('\n\033[31m\033[1m\033[3m------------ Start Jacobi method -------------')
        print('\n')

    if Zeros_Matrix is None:
        Zeros_Matrix = zeros(len(Matrix_A))

    D = diag(Matrix_A)
    R = Matrix_A - diagflat(D)

    for i in range(iteration):
        Zeros_Matrix = (Matrix_B - dot(R, Zeros_Matrix)) / D
    print(Zeros_Matrix)
    print('\n\033[31m\033[1m\033[3m------------- End Jacobi method --------------')
    return Zeros_Matrix


def gauss_saidel(Matrix_A, Matrix_B, iteration=2, accuracy=1e-6):
    n = len(Matrix_A)

    btol = linalg.norm(Matrix_B) * accuracy

    x0 = zeros(n)
    k = 0
    isActive = False
    x1 = empty(n)

    while not (isActive) and k < iteration:
        print(f'Начало итерации при k = {k}')
        print('------------------------------------')
        x1 = zeros(n)
        for i in range(n):
            x1[i] = (Matrix_B[i] - dot(Matrix_A[i, 0:i], x1[0:i]) - dot(Matrix_A[i, i + 1:n], x0[i + 1:n])) / Matrix_A[i, i]
            print("x1  = ", x1)

        r = Matrix_B - dot(Matrix_A, x1)
        isActive = (linalg.norm(r) < btol) and (linalg.norm(x1-x0) < accuracy)
        print(f'\n')
        print("x0  = ", x0)
        print("btol = %e;\nla.norm(r) = %e;\ntol = %e;\nla.norm(x1-x0) = %e;\nisActive = %s " % (btol, linalg.norm(r), accuracy, linalg.norm(x1-x0), isActive))
        x0 = x1
        print("x0  = ", x0, end='\n')
        print('------------------------------------')
        print("Окончание итерации \n\n")
        k = k + 1

    if not (isActive):
        print(f'Не сходится за {iteration} итерации.')

    return x1
