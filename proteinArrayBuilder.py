import os
import math

directory = "MFS_1"


def build_protein_arrays():

    protein_array = []
    protein_names = []
    shortest_string = math.inf

    for file in os.listdir(directory):
        f = open(os.path.join(directory, file), "r")
        protein_names.append(file)
        protein = f.read()
        protein = protein.split('\n')[1::]
        protein = ''.join(protein)
        protein_array.append(protein)

        if len(protein) < shortest_string:
            shortest_string = len(protein)

    for protein in protein_array:
        protein = protein[:shortest_string]

    return protein_array, protein_names

