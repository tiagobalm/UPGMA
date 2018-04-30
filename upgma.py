from ete3 import Tree


def abc(start, end):
    abc = []
    for i in range(ord(start), ord(end) + 1):
        abc.append(chr(i))
    return abc


def coordinates_lowest_value_matrix(matrix):
    x, y = -1, -1
    min_value = float("inf")
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] < min_value:
                min_value = matrix[i][j]
                x, y = i, j
    return x, y


def clustering_abc(abc, x, y):
    if y < x:
        x, y = y, x
    abc[x] = "(" + abc[x] + "," + abc[y] + ")"
    del abc[y]
    return abc


def update_matrix(matrix, x, y):
    if y < x:
        x, y = y, x
    row = []
    for i in range(0, x):
        row.append((matrix[x][i] + matrix[y][i]) / 2)
    matrix[x] = row

    for i in range(x + 1, y):
        matrix[i][x] = (matrix[i][x] + matrix[y][i]) / 2

    for i in range(y + 1, len(matrix)):
        matrix[i][x] = (matrix[i][x] + matrix[i][y]) / 2
        del matrix[i][y]
    del matrix[y]
    return matrix


def upgma(distance_matrix, abc):
    while len(abc) > 1:
        x, y = coordinates_lowest_value_matrix(distance_matrix)
        update_matrix(distance_matrix, x, y)
        clustering_abc(abc, x, y)
    return abc[0]


def constructing_final_tree(distance_matrix, abc):
    v = str(upgma(distance_matrix, abc)) + ";"
    t = Tree(v)
    print(t)


distance_matrix = [
    [],
    [9],
    [2, 9],
    [4, 6, 5],
    [9, 2, 9, 6],
    [10, 10, 10, 10, 10]
]

abc = abc("A", "F")

print(constructing_final_tree(distance_matrix, abc))

