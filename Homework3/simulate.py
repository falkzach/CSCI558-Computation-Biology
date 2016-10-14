#!/usr/bin/python
#!/usr/local/bin/python3

# note that the shebang is set to run from python2 instead of python3, matplotlib missing on target env for python3

from __future__ import print_function

import random
import matplotlib.pyplot as plt


FASTA_SEQUENCE_START_REGEX = '>(?P<sequence_id>.*)\n'

A = 'A'
C = 'C'
G = 'G'
T = 'T'
ACTG = [A, C, G, T]

SCORE_MATRIX_INDEX = {
    'A': 0,
    'T': 1,
    'C': 2,
    'G': 3
}
SCORE_MATRIX = [
    [4, -5, -5, -1],
    [-5, 4, -1, -5],
    [-5, -1, 4, -5],
    [-1, -5, -5, 4]
]
GAP_PENTALTY = -5

NUMBER_SIMULATIONS = 100000
SEQUENCE_LENGTH = 100
CHARACTER_PROBABILITY = 0.25


def local_score(sequence1, sequence2):
    width, height = len(sequence1) + 1, len(sequence2) + 1
    matrix = [[0 for x in range(width)] for y in range(height)]
    for x in range(1, width):
        for y in range(1, height):
            scores = [0]
            dx, dy = x - 1, y - 1
            character_x, character_y = sequence1[x - 1], sequence2[y - 1]
            ix, iy = SCORE_MATRIX_INDEX[character_x], SCORE_MATRIX_INDEX[character_y]
            local_score = SCORE_MATRIX[iy][ix]
            scores.append(matrix[dy][dx] + local_score)
            scores.append(matrix[y][dx] + GAP_PENTALTY)
            scores.append(matrix[dy][x] + GAP_PENTALTY)
            matrix[y][x] = max(scores)
    return matrix, max(map(max, matrix))


def print_matrix(matrix, X, Y):
    row_format = '{:^3}' * (len(X) + 2)
    print(row_format.format(*list('  ' + X)))
    Y = ' ' + Y
    for col, row in enumerate(matrix):
        row.insert(0, Y[col])
        print(row_format.format(*row))


if __name__ == '__main__':
    scores = []
    for x in range(NUMBER_SIMULATIONS):
        print(str(x) + '/' + str(NUMBER_SIMULATIONS), end='\r')
        X = ''.join(random.choice(ACTG) for i in range(SEQUENCE_LENGTH))
        Y = ''.join(random.choice(ACTG) for i in range(SEQUENCE_LENGTH))
        matrix, score = local_score(X, Y)
        # print_matrix(matrix, X, Y)
        # print("The score is " + str(score))
        scores.append(score)
    print(scores)

    plt.hist(scores)
    plt.title('Local Alignment Scores for ' + str(NUMBER_SIMULATIONS) + ' simulations of pairwise ' + str(SEQUENCE_LENGTH) + ' character sequences')
    plt.xlabel('Scores')
    plt.ylabel('frequency')
    plt.savefig('simulated.png')
