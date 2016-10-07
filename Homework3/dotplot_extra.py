# !/usr/local/bin/python3

import argparse
import re

FASTA_SEQUENCE_START_REGEX = '>(?P<sequence_id>.*)\n'


def get_input_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Input File")
    parser.add_argument("k_value", help="k Consecutive Matches")
    args = parser.parse_args()
    if args.input_file is None:
        print('Useage: {} <input_file> <k_value>'.format(__file__))
        exit(1)
    return args.input_file, int(args.k_value)


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
    if (match):
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


def dotplot(sequence1, sequence2, k):
    width, height = len(sequence1), len(sequence2)

    # instantiate a matrix of no matches, this allows us to only write matches
    matrix = [[' ' for x in range(width)] for y in range(height)]
    # normal scan through sequences
    for x in range(width):
        for y in range(height):
            if sequence1[x] is sequence2[y]:
                # handle edge cases
                # only the lower right segment can have hits unless k is 1
                # this optimization allows us to simply ignore everything outside of that
                match = True if (x >= k - 1 and y >= k - 1) or (k <= 1) else False

                # look back for k values to ensure it is still a match, only if we are looking for continuous matches
                if match and k > 1:
                    # check only the diagonal, nearest first
                    for i in range(1,k):
                        dx, dy = x - i, y - i
                        if dx > 0 and dy > 0:
                            if sequence1[dx] is not sequence2[dy]:
                                match = False
                                break
                if match:
                    matrix[y][x] = '.'
    return [''.join(row) for row in matrix]


def print_dotplot(sequence1, sequence2, plot):
    print(' ' + sequence1)
    for y, row in enumerate(plot):
        print(sequence2[y] + row)

if __name__ == '__main__':
    input_file, k_value = get_input_arguments()
    sequences = read_sequences_from_file(input_file)
    X = sequences['seq1']
    Y = sequences['seq2']
    plot = dotplot(X, Y, k_value)
    print_dotplot(X, Y, plot)
