import sys
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":

    match = 1
    mismatch = -4

    matches = 0
    score = 0
    
    ranges_in = sys.argv[1]
    results_in = sys.argv[2]

    ranges_df = pd.read_csv(ranges_in, header=None, delim_whitespace=True)
    results_df = pd.read_csv(results_in, header=None, delim_whitespace=True)

    found = {}
    locations = {}

    last = ranges_df.values.max()
    for i in range(0, last):
        locations[i] = -1

    for index, row in ranges_df.iterrows():
        found[index] = False
        for i in range(row[0], row[1]):
            locations[i] = index

    for index, row in results_df.iterrows():
        overlap = False
        for i in range(row[0], row[1]):
            if i < last and locations[i] != -1:
                if found[locations[i]] == False:
                    score += match
                    matches += 1
                found[locations[i]] = True
                overlap = True
                break
        if overlap == False:
            score+= mismatch

    # this is a really bad version of the scoring algorith that will run forever on larege inputs
    # for index, row in ranges_df.iterrows():
    #     ranges.append( range(row[0], row[1]) )
    #     found.append(False)

    # for index, row in results_df.iterrows():
    #     results.append( range(row[0], row[1]) )

    # for result_range in results:
    #     matched = False
    #     for (i, rna_range) in enumerate(ranges):
    #         if (len(list(set(rna_range) & set(result_range))) > 0):
    #             matched = True
    #             break

    #     if matched == True:
    #         if found[i] == True:
    #             continue
    #         found[i] = True
    #         score += match
    #         matches += 1
    #     else:
    #         score += mismatch

    print("matches: {}".format(matches))
    print("score: {}".format(score))
