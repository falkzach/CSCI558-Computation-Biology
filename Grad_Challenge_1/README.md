# Fast Fragment Assembly Program

This sample was my submission for a Graduate Programming Challenge in a Computational Biology course. The objective was to write an algorithm to assemble a series of fragmented DNA into the original sequence. It was the only submission to consistently reconstruct large sequences from real genomes and did so with a significant performance advantage. At its core it is a greedy graph problem, implemented in C++11 on a std::map and leverages OMP for parallelization.

The solution generates a graph of directed edges from a node who's suffix is the prefix of the destination node
with a weight of the number of overlapping characters.
In addition, a metagraph is kept tracking overlap weights to edges.
Those edges are then greedily consumed to construct a Hamiltonian path.
Finally, the path is walked and a result sequence constructed based on the path and tracked overlaps.
The input size is reduced by eliminating fragments that are pure substrings of other fragments.

This implementation leans on a few key assumptions

1. that the fragments provide full coverage of the original sequence

2. that fragments which are pure substrings of other fragments provide no extra information

    - while this may not be true, leaning on significant over coverage (see 1) makes this acceptable

    - working under this assumption, the size of the input set can be significantly reduced before construction of the graph

3. that greedy consumption of edges will provide an ideal reconstruction of the sequence

```
GATTA
 ATTACCAATT
    ACCAATTAC
      CAATTACC
       AATTACCAGG
          TACCAGGA
GATTACCAATTACCAGGA
```

## included files
```assemble.cpp``` - my code sample.

```fragments.txt``` - a list of short fragments produced from the sequence ```GATTACCAATTACCAGGA```.

```Makefile``` - compilation directives.

```README.md``` - this README file.

## Compile
```make```

# Run
```make check```

or

```./assemble.o ./fragments.txt```
