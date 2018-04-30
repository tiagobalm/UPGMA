from alignment import needleman_wunsch
from BLOSUM62 import *

threshold = 0


def calculate_score(alignment):
    score = 0

    for i in range(0, len(alignment[0])):
        if alignment[0][i] == alignment[2][i]:
            score += 1
        elif alignment[0][i] != '-' and alignment[2][i] != '-' and alignment[0][i] != alignment[2][i]:
            if blossum62[letters_map_proteins[alignment[0][i]]][letters_map_proteins[alignment[2][i]]] > threshold:
                score += 1

    return score


def calculate_best_score(all_alignments):
    best_score = 0

    for i in range(0, len(all_alignments)):
        score = calculate_score(all_alignments[i])
        if score > best_score:
            best_score = score

    return best_score


def build_distance_matrix(proteins):
    distance_matrix = []

    for i in range(0, len(proteins)):
        distances = []
        for j in range(0, i):
            distances.append(calculate_best_score(needleman_wunsch(proteins[i], proteins[j], -3, False)))
        distance_matrix.append(distances)

    return distance_matrix
