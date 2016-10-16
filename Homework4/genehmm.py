#!/usr/local/bin/python3

import argparse
import sys
import random
import math

A = 'A'
C = 'C'
G = 'G'
T = 'T'
ACGT = [A, C, G, T]

EMIT = 'emit'
LABEL = 'label'
MODES = [EMIT, LABEL]

START = 'start'
END = 'end'
EXON = 'exon'
INTRON = 'intron'
BACKGROUND = 'background'

P_state = {
    START: {
        EXON: 1
    },
    EXON: {
        EXON: 0.95,
        INTRON: 0.04,
        END: 0.01
    },
    INTRON: {
        INTRON: 0.95,
        EXON: 0.05
    }
}

P_emit = {
    EXON: {
        A: 0.1,
        C: 0.4,
        G: 0.4,
        T: 0.1
    },
    INTRON: {
        A: 0.4,
        C: 0.1,
        G: 0.1,
        T: 0.4
    }
}

P_BG_state = {
    START: {
        BACKGROUND: 1
    },
    BACKGROUND: {
        BACKGROUND: 0.99,
        END: 0.01
    }
}

P_BG_emit = {
    BACKGROUND: {
        A: 0.25,
        C: 0.25,
        G: 0.25,
        T: 0.25
    }
}

VITIBRI_ID = {
    0: EXON,
    1: INTRON
}

VITIBRI_LABELS = {
    EXON: 'E',
    INTRON: 'I'
}


def get_input_arguments():
    parser = argparse.ArgumentParser()
    mode_help = '|'.join(MODES)
    parser.add_argument('mode', help='HMM mode, ' + mode_help)
    args = parser.parse_args()
    if args.mode is None or args.mode not in MODES:
        print('Useage: {} [mode <{}>] '.format(__file__,  mode_help))
        exit(1)
    return args.mode


# http://stackoverflow.com/a/3679747
def weighted_choice(choices):
    total = sum(w for c, w in choices)
    r = random.uniform(0, total)
    upto = 0
    for c, w in choices:
        if upto + w >= r:
            return c
        upto += w
    assert False, "Shouldn't get here"


def change_state(current_state):
    P = P_state[current_state].items()
    return weighted_choice(P)


def emit_character(current_state):
    P = P_emit[current_state].items()
    char = weighted_choice(P)
    sys.stdout.write(char)


def emit_sequence():
    state = START
    while True:
        state = change_state(state)
        if state == END:
            break
        emit_character(state)
    sys.stdout.write('\n')


def read_stdin():
    return sys.stdin.read().rstrip()


def label_sequence(data):
    X = [x for x in data]
    label = ''
    width, height = len(X), 2
    matrix = [[0 for x in range(width)] for y in range(height)]


    character = X[0]
    matrix[0][0] = 1 * P_emit[EXON][character]# start position, this should be multiplied by first char, change, state, emit, remember
    # is matrix[1][0] 0 since it cant start?
    for x in range(1, width):
        for y in range(height):
            state_key = VITIBRI_ID[y]
            P_s = P_state[state_key]
            P_e = P_emit[state_key]
            # P(sequence | gene model)/P(sequence | background model)
            # log2(p) FOR LOG ODDS RATIO
            character = X[x]
            # p = math.log(
            #     (P_s[state_key] * P_e[character]) /
            #     (P_BG_state[BACKGROUND][BACKGROUND] * P_BG_emit[BACKGROUND][character]),
            #     2
            # )
            p = P_s[state_key] * P_e[character]
            e = matrix[0][x-1] * p
            i = matrix[1][x-1] * p
            v = max(e, i)
            matrix[y][x] = v

        # pretty sure this is very wrong
        winning_state_key = EXON if matrix[0][x] > matrix[1][x] else INTRON
        label += VITIBRI_LABELS[winning_state_key]

    print(data)
    print(label)
    for row in matrix:
        print(row)
    print(matrix[0][width - 1])
    print(matrix[1][width - 1])


if __name__ == '__main__':
    mode = get_input_arguments()
    if mode == EMIT:
        emit_sequence()
    elif mode == LABEL:
        data = read_stdin()
        label_sequence(data)
    exit(0)
