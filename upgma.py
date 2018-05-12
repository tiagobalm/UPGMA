from ete2 import Tree, TreeStyle
from matrixBuilder import build_distance_matrix
from proteinArrayBuilder import build_protein_arrays
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


def coordinates_lowest_value_matrix(matrix):
    x, y = -1, -1
    min_value = float("inf")
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] < min_value:
                min_value = matrix[i][j]
                x, y = i, j
    return x, y


def clustering_abc(matrix, abc, x, y):
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
        clustering_abc(distance_matrix, protein_labels, x, y)
        if len(distance_matrix) == 2:
            distance_to_root = distance_matrix[1][0] / 2
    return protein_labels[0]


def constructing_final_tree(distance_matrix, protein_labels):
    v = str(upgma(distance_matrix, protein_labels)) + ";"
    t = Tree(v)
    ts = TreeStyle()
    ts.show_leaf_name = True
    t.convert_to_ultrametric()
    ts.scale = 120
    t.show(tree_style=ts)


proteinsArray, proteinNames = build_protein_arrays()

logging.debug("Building distance matrix. Please wait...")

distanceMatrix = build_distance_matrix(proteinsArray)

logging.debug("Building tree....")

print(constructing_final_tree(distanceMatrix, proteinNames))

