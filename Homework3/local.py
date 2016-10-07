# !/usr/local/bin/python3

import argparse
import re

FASTA_SEQUENCE_START_REGEX = '>(?P<sequence_id>.*)\n'

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


def get_input_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Input File")
    args = parser.parse_args()
    if args.input_file is None:
        print('Useage: {} <input_file>'.format(__file__))
        exit(1)
    return args.input_file


def get_match(line):
    pattern = re.compile(FASTA_SEQUENCE_START_REGEX)
    return pattern.match(line)


def is_sequence_start(line):
    match = get_match(line)
    if match:
        return True
    return False


def get_sequence_number(line):
    match = get_match(line)
    if match:
        return match.group('sequence_id')


def add_line_to_sequence_dictionary(sequence_id, line, sequences):
    line = line.rstrip()
    if sequence_id not in sequences:
        sequences[sequence_id] = line
    else:
        sequences[sequence_id] += line


def read_sequences_from_file(file_name):
    sequence_id = None
    sequences = {}

    f = open(file_name, 'r')
    for line in f:
        if is_sequence_start(line):
            sequence_id = get_sequence_number(line)
        elif sequence_id is not None:
            add_line_to_sequence_dictionary(sequence_id, line, sequences)
    return sequences


def score_dynamic_programming(sequence1, sequence2):
    width, height = len(sequence1) + 1, len(sequence2) + 1
    matrix = [[0 for x in range(width)] for y in range(height)]
    for x in range(1, width):
        for y in range(1, height):
            scores = []
            dx, dy = x - 1, y - 1
            ix, iy = SCORE_MATRIX_INDEX[sequence1[x - 1]], SCORE_MATRIX_INDEX[sequence2[y - 1]]
            local_score = SCORE_MATRIX[iy][ix]
            scores.append(matrix[dy][dx] + local_score)
            scores.append(matrix[y][dx] + GAP_PENTALTY)
            scores.append(matrix[dy][x] + GAP_PENTALTY)
            matrix[y][x] = max(scores)
    return matrix, matrix[height - 1][width - 1]


if __name__ == '__main__':
    input_file = get_input_arguments()
    sequences = read_sequences_from_file(input_file)
    X = sequences['seq1']
    Y = sequences['seq2']

    matrix, score = score_dynamic_programming(X, Y)
    print(score)
    print(X)
    print(Y)
    for row in matrix:
        print(row)
