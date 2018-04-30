from ete3 import Tree
from matrixBuilder import build_distance_matrix

protein1 = "ARNDCQEGHILKMFPSTWYVBZX"
protein2 = "ARHDCAEGHILKMFPSTWYVBZX"
protein3 = "ARNDCAEGHCLKMFPSTWYVBZX"
protein4 = "ARNDCQEGHCLKMFPSTWYVBZX"
protein5 = "ARNDCQEGHILKMFPSTWYWBZX"
protein6 = "ARNDCQEGHILKMFPSTWYVWWW"

proteinsArray = [protein1, protein2, protein3, protein4, protein5, protein6]
proteinNames = ["protein1", "protein2", "protein3", "protein4", "protein5", "protein6"]


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


def upgma(distance_matrix, protein_labels):
    while len(protein_labels) > 1:
        x, y = coordinates_lowest_value_matrix(distance_matrix)
        update_matrix(distance_matrix, x, y)
        clustering_abc(protein_labels, x, y)
    return protein_labels[0]


def constructing_final_tree(distance_matrix, protein_labels):
    v = str(upgma(distance_matrix, protein_labels)) + ";"
    t = Tree(v)
    return t


distanceMatrix = build_distance_matrix(proteinsArray)
print(constructing_final_tree(distanceMatrix, proteinNames))
