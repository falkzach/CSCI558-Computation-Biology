#!/usr/local/bin/python3

import argparse
import subprocess
import re

def get_input_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('sto_file', help='Input File')
    parser.add_argument('fasta_db', help='Output File')
    parser.add_argument('threshold', help='E-value Threshold')
    args = parser.parse_args()
    if args.sto_file is None or args.fasta_db is None or args.threshold is None:
        print('Useage: {} [sto_file] [fasta_db] '.format(__file__))
        exit(1)
    return args.sto_file, args.fasta_db, args.threshold


if __name__ == '__main__':
    file = 'hmmdata'
    sto_file, fasta_db, threshold = get_input_arguments()
    call1 = ['hmmbuild', file, sto_file]
    process1 = subprocess.Popen(call1, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process1.communicate()

    call2 = ['hmmsearch', '-E {}'.format(threshold), '--noali', file, fasta_db]
    process2 = subprocess.Popen(call2, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process2.communicate()
    clean = [line.decode('UTF-8') for line in out.splitlines()]

    # nasty poorly contrived regexes
    query_regex = '>>\s+(?P<sequence>.*)'
    value_regex = '\s+\d\s+!\s+(?P<score>\S+)\s+(?P<bias>\S+)\s+.*'
    query_pattern = re.compile(query_regex)
    value_pattern = re.compile(value_regex)

    count = 0
    match = False
    sequence = ''
    results = []
    # this is diry, so have some comments to explain
    # parse out data from output
    for line in clean:
        # if we've found a sequence line, we need to skip some lines
        if count >= 0:
            count -= 1
            if count != 0:
                continue
        # otherwise check for a new sequence string
        else:
            match = query_pattern.match(line)

        # if we've found a sequence line and then skipped formatting lines, look at scores
        if match and (count == 0):
            values = value_pattern.match(line)
            # extract eval and bias from the score line via regex
            score = values.group('score')
            bias = values.group('bias')
            results.append((sequence, score, bias))
            match = False

        # if we've found a line, set the counter to skip formatting lines
        if match:
            sequence = match.group('sequence')
            count = 3

    for result in results:
        print('{}\t{}   {}'.format(result[0], result[1], result[2]))

