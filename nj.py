from ete2 import Tree, TreeStyle, random_color
from matrixBuilder import build_distance_matrix
from proteinArrayBuilder import build_protein_arrays
import logging

       
def symmetric_matrix(matrix):
        dimension=len(matrix)
        symmetric_matrix=[[0 for x in range(dimension)] for y in range(dimension)]
        for i in range(dimension):
                for j in range(dimension):
                    if i!=j and i<j:
                        symmetric_matrix[i][j]=matrix[j][i]
                    if i!=j and i>j:
                        symmetric_matrix[i][j]=matrix[i][j]
        return symmetric_matrix
            
       
def r(matrix):
    l=len(matrix)
    r=[0]*l
    for i in range(l):
        r[i]=sum(matrix[i])/(l-2)
    return r


def D_i_j(matrix):
    vector=r(symmetric_matrix(matrix))
    D=matrix
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i!=j and i>j:
                D[i][j]=matrix[i][j]-vector[i]-vector[j]
    return D

                    
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
        row.append((matrix[x][i] + matrix[y][i]-matrix[y][x]) / 2)
    matrix[x] = row

    for i in range(x + 1, y):
        matrix[i][x] = (matrix[i][x] + matrix[y][i]-matrix[y][x]) / 2

    for i in range(y + 1, len(matrix)):
        matrix[i][x] = (matrix[i][x] + matrix[i][y]-matrix[y][x]) / 2
        del matrix[i][y]
    del matrix[y]
    return matrix




def neighbor_joining(distance_matrix, protein_labels):
    while len(distance_matrix) > 2:
       x, y = coordinates_lowest_value_matrix(D_i_j(distance_matrix))
       update_matrix(distance_matrix, x, y)
       clustering_abc(protein_labels, x, y)
    if len(distance_matrix)==2:
            clustering_abc(protein_labels, 0, 1)
    return protein_labels[0]



def my_layout(node):
    if node.is_leaf():
       node.img_style["size"] = 4
       node.img_style["shape"] = "sphere" 
       node.img_style["fgcolor"] = "black"
    else:
       node.img_style["size"] = 0

       
def constructing_final_tree(distance_matrix, protein_labels):
    v = str(neighbor_joining(distance_matrix, protein_labels)) + ";"
    t = Tree(v)
    t.dist = 0
    ts = TreeStyle()
    ts.mode = "c" 
    ts.show_leaf_name = True
    ts.layout_fn = my_layout
    t.show(tree_style=ts)


proteinsArray, proteinNames = build_protein_arrays()

logging.debug("Building distance matrix. Please wait...")

distanceMatrix = build_distance_matrix(proteinsArray)

logging.debug("Building tree....")


print(constructing_final_tree(distanceMatrix, proteinNames))






    

    
        
         
    

    
    
