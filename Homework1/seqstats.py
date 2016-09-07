#!/usr/local/bin/python3

import argparse
import re
import copy
import collections

FASTA_SEQUENCE_START_REGEX ='>(?P<sequence_id>.*)\n'
A = 'A'
C = 'C'
G = 'G'
T = 'T'
ACTG = [A,C,G,T]
LEN = 'Len'
COUNTS = {LEN: 0, A: 0, C: 0, T: 0, G: 0}
TABLE_FORMAT = '{:^9} {:^3} {:^3} {:^3} {:^3} {:^3}'

def getInputArguments():
    parser = argparse.ArgumentParser()
    #-i <inputFile>
    parser.add_argument("-i", "--inputFile", help="Input File")
    args = parser.parse_args()
    if args.inputFile is None:
        print('Useage: {} -i <inputFile>'.format(__file__))
        exit(1)
    return args.inputFile

def getMatch(line):
    pattern = re.compile(FASTA_SEQUENCE_START_REGEX)
    return pattern.match(line)

def isSequenceStart(line):
    match = getMatch(line)
    if(match):
        return True
    return False

def getSequenceNumber(line):
    match = getMatch(line)
    if(match):
        return match.group('sequence_id')

def addLineToSequenceToDictionary(sequence_id, line, sequences):
    line = line.rstrip()
    if(sequence_id not in sequences):
        sequences[sequence_id] = line
    else:
        sequences[sequence_id] += line

def readFile(inputFile):
    sequence_id = None
    sequences = {}

    f = open(inputFile, 'r')
    for line in f:
        if isSequenceStart(line):
            sequence_id = getSequenceNumber(line)
        elif(sequence_id is not None):
            addLineToSequenceToDictionary(sequence_id, line, sequences)
    return sequences

def countSequences(sequences):
    return len(sequences)

def countOccurancesinSequence(sequence):
    counts = copy.deepcopy(COUNTS)
    for character in sequence:
        if character in ACTG:
            counts[character] += 1
            counts[LEN] += 1
    return counts

def countOccurancesinSequences(sequences):
    sequenceCounts = {}
    for sequence_id, sequence in sequences.items():
        sequenceCounts[sequence_id] = countOccurancesinSequence(sequence)
    return sequenceCounts

def analizeSequences(sequences):
    return countSequences(sequences), countOccurancesinSequences(sequences)

if __name__ == '__main__':

    inputFile = getInputArguments()

    sequences = readFile(inputFile)

    numberOfSequences, sequenceCounts = analizeSequences(sequences)

    orderedSequenceCounts = collections.OrderedDict(
        sorted(sequenceCounts.items())
    )
    print('The file contains ' + str(numberOfSequences) + ' sequences.')
    print(TABLE_FORMAT.format('Seq_name',LEN,A,C,G,T))
    print(TABLE_FORMAT.format('--------','---','---','---','---','---'))
    for sequence_id, sequenceCounts in orderedSequenceCounts.items():
        sequenceName = str(sequence_id)
        print(TABLE_FORMAT.format(
            sequenceName,
            sequenceCounts[LEN],
            sequenceCounts[A],
            sequenceCounts[C],
            sequenceCounts[G],
            sequenceCounts[T]
        ))
