# !/usr/local/bin/python3

import argparse
import re

FASTA_SEQUENCE_START_REGEX = '>(?P<sequence_id>.*)\n'


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


def score_alignment(sequence1, sequence2):
    score = 0
    for charactere1, charactere2 in zip(sequence1, sequence2):
        if charactere1 is charactere2:
            score += 1
        else:
            score -= 1
    return score


if __name__ == '__main__':
    input_file = get_input_arguments()
    sequences = read_sequences_from_file(input_file)
    seq1 = sequences['seq1']
    seq2 = sequences['seq2']
    score = score_alignment(seq1, seq2)
    print(score)

