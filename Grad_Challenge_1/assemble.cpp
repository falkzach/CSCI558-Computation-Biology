/**
 * Zachary Falkner
 * CSCI558 - Computational Biology
 * Assembly
 * 
 * My implementation for this solution generates a graph of directed edges from a node who's suffix is the prefix of the destination node
 * with a weight of the number of overlapping characters.
 * In addition a metagraph is kept tracking overlap weights to edges.
 * Those edges are then greedily consumed to construct a hamiltonian path.
 * Finally the path is walked and a result sequence constructed based on the path and tracked overlaps.
 * The input size is reduce by eliminating fragments that are pure substrigns of other fragments.
 * 
 * Open MP pragmas are included to signifigantly increase performance on large input sets.
 * 
 * This implementation leans on a few key assumptions
 * 1.) that the fragments provide full coverage of the original sequence
 * 2.) that fragments which are pure substrings of other fragments provide no extra information
 *      - while this may not be true, leaning on signifigant over coverage (see 1) makes this acceptable
 *      - working under this assumption, the size of the input set can be signfigently reduce before construction of the graph
 * 3.) that greedy consumption of edges will provide an ideal reconstruction of the sequence
 * 
 * It does alright.
 * 
 * GATTA
 *  ATTACCAATT
 *     ACCAATTAC
 *       CAATTACC
 *        AATTACCAGG
 *           TACCAGGA
 * GATTACCAATTACCAGGA
 * 
**/

/**
 * GRAD CHALLENGE 1 (GENOME ASSEMBLY): Create a program to try to assemble a linear chromosome (a string of 'G', 'A', 'T', and 'C' [uppercase only] using only short reads). Your program must be entirely original and must not use any external libraries (except for numpy in python or the C++ standard library in C++). Your program must consist of a single source file (if using python) or compile from a single source file without any extra compilation flags (if using C++). C++ programs will be compiled with -std=c++11 (and -O3 and -march=native and -fopenmp), so using the C++11 standard is fine. Regardless of language, your program should accept one command line parameter: the name of a text file containing all of the short reads (one short read per line). Your program should output only one thing: the string ('G', 'A', 'T', and 'C' [uppercase only]) that you believe was used to generate the reads (whitespace in your output does not matter); this is your genome assembly. Your submission must run in <10 minutes on <10000 fragments each with length >10bp and <100bp.
 * 
 * Here is the program that will be used to generate short reads: make_fragments.py
 * To learn how to use make_fragments.py, run it without arguments.
 * Calling
 * 
 * python make_fragments.py 'GATTACCAATTACCAGGA' 20 5 10
 * resulted in output
 * TTACCA
 * CCAGGA
 * ATTACCAATT
 * TACCAGGA
 * GATTA
 * GATTA
 * CCAGG
 * ATTAC
 * CCAGGA
 * CAATTACC
 * AATTACCAG
 * ACCAGGA
 * TACCAGG
 * ACCAATTAC
 * CCAATT
 * CCAATTA
 * CCAGGA
 * ACCAGGA
 * TACCAG
 * AATTACCAGG
 * Your program is supposed to read a file of the above output and try to infer the original string, ``GATTACCAATTACCAGGA''.
 * 
 * The student with the best alignment to my secret, original string will be decided using Needleman-Wunsch with parameters: gap=-4, mismatch=-2, match=1. Y.
 * 
 * If you want to improve your performance when testing, here is a C++ Needleman-Wunsch implementation and a python program to draw a heatmap from your pass-through score matrix (which the C++ program outputs to stdout).
 * 
 * DUE DATE: October 8, 2018
**/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

#define DISPLAY_ALIGNMENTS 0
#define MIN_FRAGMENT_OVERLAP 3 //TODO: BUG: a min overlap of 2 leads to a core dump.

/* graphLike is a map of node id's to vectors of pairs of node id and weight */
using graphLike = std::map<int, std::vector<std::pair<int, int>>>;

/* read each line from the source file in as a fragment, cleaning newlines and returns */
void read_fragments(std::vector<std::string> & fragments, const std::string & path) {
        std::string fragment;
        std::ifstream fragments_fin(path);
        while(std::getline(fragments_fin, fragment)) {
            /* clean newlines and returns */
            fragment.erase(std::remove(fragment.begin(), fragment.end(), '\n'), fragment.end());
            fragment.erase(std::remove(fragment.begin(), fragment.end(), '\r'), fragment.end());
            fragments.push_back(fragment);
        }
        fragments_fin.close();
}

/* remove all fragments that are substrings of other fragments */
void prune_substring_fragments(std::vector<std::string> & fragments) {
    std::set<int> prune_fragments;
    unsigned int I, J;
    #pragma omp parallel for private(J) schedule(dynamic, 3)
    for(I=0; I < fragments.size(); ++I) {
        for (J=I+1; J < fragments.size(); ++J) {
            if(fragments[I].length() <= fragments[J].length()) {
                /* check if I is substring of J, erase I if so */
                if(fragments[J].find(fragments[I]) != std::string::npos) {
                    #pragma omp critical
                    prune_fragments.insert(I);
                }
            } else {
                /* check if J is substring of I, erase J if so */
                if(fragments[I].find(fragments[J]) != std::string::npos) {
                    #pragma omp critical
                    prune_fragments.insert(J);
                }
            }
        }
    }
    /* remove fragments from greatest index to smallest, don't have to change other indexes to prune as vector shrinks */
    for(std::set<int>::reverse_iterator rit=prune_fragments.rbegin(); rit != prune_fragments.rend(); rit++) {
        int id = *rit;
        fragments.erase(fragments.begin() + id);
    }
}

/* add an edge to the scoreGraph that points to the edge in the graph */
void add_score_edge(graphLike & scoreGraph, const int score, const int & A, const int & B) {
    #pragma omp critical
    if (scoreGraph.count(score) == 0) {
        std::vector<std::pair<int, int>> edges;
        scoreGraph[score] = edges;
    }
    std::pair<int, int> edge = {A, B};
    #pragma omp critical
    scoreGraph[score].push_back(edge);
}

/**
 * build out a directed graph where an edge represents an overlap between the 2 strings
 * that is, the source node has a sufix that is a prefix to the destination node with some overlap
**/
void build_graph(const std::vector<std::string> fragments, graphLike & graph, graphLike & scoreGraph) {
    /* put all nodes in the graph with no edges */
    #pragma omp parallel for
    for(unsigned long i=0; i < fragments.size(); ++i) {
        std::vector<std::pair<int, int>> adjacencies;
        #pragma omp critical
        graph[i] = adjacencies;
    }

    unsigned long I, J;
    /* build out edges of graph */
    #pragma omp parallel for private(J) schedule(dynamic, 3)
    for(I=0; I < fragments.size(); ++I) {
        for (J=I+1; J < fragments.size(); ++J) {
            int weightI = 0;  //suffix of I prefix of J I -> J
            int weightJ = 0;  //suffix of J prefix of I J -> I
            int bestI = 0;
            int bestJ = 0;
            for (unsigned int i = 0, j=fragments[J].length() - 1; i < fragments[I].length() && j >= 0; ++i, --j) {
                    weightJ = 0;
                    for (unsigned int c = 0, d = j; c <= i && d < fragments[J].length(); ++c, ++d) {
                        if(fragments[I][c] == fragments [J][d]) {
                            ++weightJ;
                        } else {
                            weightJ = 0;
                            break;
                        }
                    }
                    /* keep the best scoring overlap */
                    bestJ = std::max(bestJ, weightJ);
            }

            for (unsigned int i=fragments[I].length() -1, j=0; i >=0 && j < fragments[J].length(); --i, ++j) {
                    weightI = 0;
                    for(unsigned int c = i, d = 0; c < fragments[I].length() && d <= j; ++c, ++d) {
                        if(fragments[I][c] == fragments [J][d]) {
                            ++weightI;
                        } else {
                            weightI = 0;
                            break;
                        }
                    }
                    /* keep the best scoring overlap */
                    bestI = std::max(bestI, weightI);
            }

            /* add edges to graph */

            if (bestJ >= MIN_FRAGMENT_OVERLAP) {
                std::pair<int, int> edge = {I, bestJ};
                #pragma omp critical
                graph[J].push_back(edge);
                add_score_edge(scoreGraph, bestJ, J, I);
            }

            if (bestI >= MIN_FRAGMENT_OVERLAP) {
                std::pair<int, int> edge = {J, bestI};
                #pragma omp critical
                graph[I].push_back(edge);
                add_score_edge(scoreGraph, bestI, I, J);
            }

        }
    }
}

/* a utility function to output a graphs adjacency lists */
void print_adjacency_lists(const graphLike & graph){
    for (auto it=graph.begin(); it!=graph.end(); ++it) {
        std::cout << it->first << ":\t";
        for (auto pair=it->second.begin(); pair!=it->second.end(); ++pair) {
            std::cout << pair->first << "(" << pair->second << "),\t";
        }
        std::cout << std::endl;
    }
}

/* append a fragment to the front of the result with some overlap, not used*/
void append_prefix_string(std::string & result, const std::string & P, const int & overlap) {
    result = P.substr(0, overlap) + result;
}

/* append a fragment to the end of the result with some overlap*/
void append_suffix_string(std::string & result, const std::string & S, const int & overlap) {
    result += S.substr(overlap, S.length());
}

/* create a greedy Hamiltonian(?) path and walk it */
std::string assemble(const std::vector<std::string> & fragments, const graphLike & scoreGraph) {
    std::string result;
    std::map<int, int> path;
    std::map<int, int> overlaps;
    std::set<int> visited;
    int current_node;

    /* greadily add edges to a path by the highest weight */
    for (auto score_rit=scoreGraph.rbegin(); score_rit!=scoreGraph.rend(); ++score_rit) {
        int overlap = score_rit->first;
        /* for each edge with the current weight */
        for(unsigned int i=0; i < score_rit->second.size(); ++i) {
            std::pair<int, int> edge = score_rit->second[i];
            /* if the destination node has not been visited and the source node is not already a part of the path */
            if(visited.count(edge.second) == 0  && path.count(edge.first) == 0) {
                /* mark the destination node as having been visited */
                visited.insert(edge.second);
                /* and add the edge to the path */
                path[edge.first] = edge.second;
                /* track the overlap for easy lookup */
                overlaps[edge.first] = overlap;
            }
        }
    }

    /* find the first node that was not visited, this should be the start */
    for(unsigned int i=0; i< fragments.size(); ++i) {
        if(visited.count(i) == 0) {
            current_node = i;
            break;
        }
    }

    /* set the result to be the first fragment */
    result = fragments[current_node];

    int overlap, offset;
    if (DISPLAY_ALIGNMENTS == 1) {
        offset = result.length() - overlaps[current_node];
        std::cout << result << std::endl;
    }
    /* while there is a next node, lookup the overlap and append the fragment */
    while (path.count(current_node) > 0) {
        overlap = overlaps[current_node];
        /* progress the current node */
        current_node = path[current_node];
        /* append to result */
        append_suffix_string(result, fragments[current_node], overlap);
        /* this will print staggered fragments with alignment */
        if (DISPLAY_ALIGNMENTS == 1) {
            for(int o=0; o<offset;++o) { std::cout << " "; }
            std::cout << fragments[current_node] << std::endl;
            offset += fragments[current_node].length() - overlaps[current_node];
        }
    }
    return result;
}

int main(int argc, char**argv) {
    if (argc == 2) {
        std::vector<std::string> fragments;
        graphLike graph;
        graphLike scoresGraph; //std::map<score, std::vector <std::pair<to, from>>>
        std::string result;

        /* read in fragments */
        read_fragments(fragments, argv[1]);

        /* prune pure substrings as they add no new information (Assumption) */
        prune_substring_fragments(fragments);

        /* Build Graph */
        build_graph(fragments, graph, scoresGraph);

        // print_adjacency_lists(graph);
        
        /* Assemble... */
        result = assemble(fragments, scoresGraph);

        /* Output */
        std::cout << result << std::endl;
        return 0;
    }
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
}
