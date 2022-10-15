from IPython.display import display
from sympy import *

init_printing()


def find_basis(m):
    rows, cols = m.shape
    basis = [None] * rows
    for i in range(1, cols):
        temp = list(m[:, i])
        if temp.count(1) == 1 and temp.count(0) == rows - 1:
            basis[temp.index(1)] = i
    print("Базисные векторы: ", basis)
    return basis


def simplex(A, coeffs, names):
    rows, cols = A.shape
    basis = find_basis(A)
    if None in basis:
        print('Нет базиса')
        return
    bases_names = Matrix([names[i] for i in basis])
    basis_coeffs = Matrix([coeffs[i] for i in basis])
    A01 = bases_names.row_join(basis_coeffs)
    feasible_solution = [0] * cols
    for i in range(rows):
        feasible_solution[basis[i]] = A[i, 0]
    print(f"Опорное решение: {feasible_solution[1:]}")

    deltas = []
    for j in range(A.cols):
        x = coeffs[j] - A[:, j].dot(basis_coeffs)
        deltas.append(x)

    A_title = Matrix([["Базис", "C_баз"] + names])
    A_rows = A01.row_join(A)
    A_delta = Matrix([["\u3164", "\u3164"] + deltas])
    # итоговая таблица
    A_all = A_title.col_join(A_rows).col_join(A_delta)
    display(A_all)

    while any(x > 0 for x in deltas[1:]) or deltas[1:].count(0) > rows:

        # запрос продолжения для вырожденной задачи
        if deltas[1:].count(0) > rows:
            go_on = input("Вырожденный случай. Продолжить решение? (y/n)")
            if go_on != "y": break

        unlimited = false
        for j in range(1, A.cols):
            if deltas[j] > 0 and all(x <= 0 for x in A[:, j]):
                print("Целевая функция не ограничена!")
                unlimited = true
                break
        if unlimited:
            break

        j_cur = int(input("укажите номер столбца, вводимого в базис: "))
        i_cur = int(input("укажите номер строки с удаляемым из базиса вектором: "))
        A_pivot = A[i_cur, j_cur]  # ведущий элемент

        for i in range(A.rows):
            if i == i_cur:
                continue
            A.row_op(i, lambda var, j:
            (var * A_pivot - A[i_cur, j] * A[i, j_cur]) / A_pivot
                     )

        A.row_op(i_cur, lambda var, j: var / A_pivot)

        basis = find_basis(A)
        bases_names = Matrix([names[i] for i in basis])
        basis_coeffs = Matrix([coeffs[i] for i in basis])
        A01 = bases_names.row_join(basis_coeffs)

        deltas = []
        for j in range(A.shape[1]):
            x = coeffs[j] - A[:, j].dot(basis_coeffs)
            deltas.append(x)

        A_title = Matrix([["Базис", "C_баз"] + names])
        A_rows = A01.row_join(A)
        A_delta = Matrix([["\u3164", "\u3164"] + deltas])
        A_all = A_title.col_join(A_rows).col_join(A_delta)

        feasible_solution = [0] * cols
        for i in range(rows):
            feasible_solution[basis[i]] = A[i, 0]
        print(f"\nОпорное решение: {feasible_solution[1:]}")

        display(A_all)

    return (A)


def main():
    c = Matrix([[0, 15, 6, 12, 24, 0, 0, 0]])
    A = Matrix([[1200, 4, 2, 1, 8, 1, 0, 0],
                [600, 2, 10, 6, 0, 0, 1, 0],
                [1500, 3, 0, 6, 1, 0, 0, 1]])

    names = [f'A{i}' for i in range(c.shape[1])]
    simplex(A, c, names)


if __name__ == '__main__':
    main()
