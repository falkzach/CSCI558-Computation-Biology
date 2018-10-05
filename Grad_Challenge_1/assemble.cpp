/** 
 * GRAD CHALLENGE 1 (GENOME ASSEMBLY): Create a program to try to assemble a linear chromosome (a string of 'G', 'A', 'T', and 'C' [uppercase only] using only short reads). Your program must be entirely original and must not use any external libraries (except for numpy in python or the C++ standard library in C++). Your program must consist of a single source file (if using python) or compile from a single source file without any extra compilation flags (if using C++). C++ programs will be compiled with -std=c++11, so using the C++11 standard is fine. Regardless of language, your program should accept one command line parameter: the name of a text file containing all of the short reads (one short read per line). Your program should output only one thing: the string ('G', 'A', 'T', and 'C' [uppercase only]) that you believe was used to generate the reads (whitespace in your output does not matter); this is your genome assembly. Your submission must run in <10 minutes on <10000 fragments each with length >10bp and <100bp.
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
 * The student with the best alignment to my secret, original string will be decided using Needleman-Wunsch with parameters: gap=-4, mismatch=-2, match=1.
 * 
 * DUE DATE: October 8, 2018
 * 
**/

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

/* graphLike is a map of node id's to vectors of pairs of node id and weight */
using graphLike = std::map<int, std::vector<std::pair<int, int>>>;

/**
 * G-> 0 = 00
 * A-> 1 = 01
 * T-> 2 = 10
 * C-> 3 = 11
 * -1uc otherwise
**/ 
const unsigned char nuc_to_code[] = {255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 1, 255, 3, 255, 255, 255, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};
const unsigned char code_to_nuc[] = {'G', 'A', 'T', 'C'};

/**
 * provides a mask for the length of the fragment, keeps n bitpacked charters where n is the index
**/ 
const unsigned long nibs_to_mask[] = {0x0, 0xc000000000000000, 0xf000000000000000, 0xfc00000000000000, 0xff00000000000000, 0xffc0000000000000, 0xfff0000000000000, 0xfffc000000000000, 0xffff000000000000, 0xffffc00000000000, 0xfffff00000000000, 0xfffffc0000000000, 0xffffff0000000000, 0xffffffc000000000, 0xfffffff000000000, 0xfffffffc00000000, 0xffffffff00000000, 0xffffffffc0000000, 0xfffffffff0000000, 0xfffffffffc000000, 0xffffffffff000000, 0xffffffffffc00000, 0xfffffffffff00000, 0xfffffffffffc0000, 0xffffffffffff0000, 0xffffffffffffc000, 0xfffffffffffff000, 0xfffffffffffffc00, 0xffffffffffffff00, 0xffffffffffffffc0, 0xfffffffffffffff0, 0xfffffffffffffffc, 0xffffffffffffffff, };

constexpr unsigned int NUCLEOTITES_PER_BLOCK = sizeof(unsigned long) * 8 / 2;

void bit_pack(unsigned long * & result, unsigned long & num_blocks, const std::string & genome) {

	num_blocks = genome.size() / NUCLEOTITES_PER_BLOCK;

	if (num_blocks * NUCLEOTITES_PER_BLOCK != genome.size()) {
		++num_blocks;
	}

	result = (unsigned long*)calloc(num_blocks, sizeof(unsigned long));

	unsigned long i = 0;

	for (unsigned long block=0; block<num_blocks-1; ++block) {
		unsigned long next_block = 0;
		for (unsigned char j=0; j<NUCLEOTITES_PER_BLOCK; ++j, ++i) {
			next_block <<= 2;

			next_block |= nuc_to_code[ genome[i] ];
		}
		result[block] = next_block;
	}

	// this loop is the inner body of the last itteration of the above loop
	unsigned long next_block = 0;
	for (unsigned char j=0; j<NUCLEOTITES_PER_BLOCK; ++j, ++i) {
		next_block <<= 2;

		if (i < genome.size()) {
			next_block |= nuc_to_code[ genome[i] ];
		}
	}
	result[num_blocks-1] = next_block;
}

//TODO: FIXME!
std::string unpack(const unsigned long * genome, const unsigned long & num_blocks, unsigned long length) {
    std::string result;
    unsigned long chars_left = length;
    for (unsigned long block=0; block < num_blocks - 1; ++block) {
        unsigned long current_block = genome[block];
        for (unsigned char j=0; j<NUCLEOTITES_PER_BLOCK; ++j, --chars_left) {
                    unsigned long nib = current_block & nibs_to_mask[1];
        nib >>= (NUCLEOTITES_PER_BLOCK * 2) - 2;
        char c = code_to_nuc[nib];
        result += c;
        current_block <<= 2;
        }
	}

    unsigned long current_block = genome[num_blocks - 1];
    for (unsigned char j=0; j<NUCLEOTITES_PER_BLOCK; ++j, --chars_left) {
        if(chars_left <= 0)
        {
            break;
        }
        unsigned long nib = current_block & nibs_to_mask[1];
        nib >>= (NUCLEOTITES_PER_BLOCK * 2) - 2;
        char c = code_to_nuc[nib];
        result += c;
        current_block <<= 2;
    }
	return result;
}

void print_bitpacked_string(const unsigned long* genome, const unsigned long & num_blocks) {
    std::cout << "packed: ";
	for (unsigned long block=0; block<num_blocks; ++block) {
		std::cout << std::hex << genome[block];
	}
	std::cout << std::endl;
}

void add_score_edge(graphLike & scoreGraph, const int score, const int A, const int B) {
    if (scoreGraph.count(score) == 0) {
        std::vector<std::pair<int, int>> edges;
        scoreGraph[score] = edges;
    }
    std::pair<int, int> edge = {A, B};
    scoreGraph[score].push_back(edge);
    // std::cout << "New Score Edge: " << A << " -> " << B << " score: " << score << std::endl;
}

/* String version TODO: bitpack*/
void build_graph(const std::vector<std::string> fragments, graphLike & graph, graphLike & scoreGraph) {
    /* put all nodes in the graph with no edges */
    for(unsigned long i=0; i < fragments.size(); ++i) {
        std::vector<std::pair<int, int>> adjacencies;
        graph[i] = adjacencies;
    }

    /* build out edges of graph */
    unsigned long I, J;
    // #pragma omp parallel for private(J)
    for(I=0; I < fragments.size(); ++I) {
        for (J=I+1; J < fragments.size(); ++J) {
                int weightI = 0;  //suffix of I prefix of J I -> J
                int weightJ = 0;  //suffix of J prefix of I J -> I

                for (int i = 0, j=fragments[J].length() - 1;
                    i < fragments[I].length() && j >= 0;
                    ++i, --j) {
                        weightJ = 0;
                        for (int c = 0, d = j; c < i && d < fragments[J].length(); ++c, ++d) {
                            if(fragments[I][c] == fragments [J][d]) {
                                ++weightJ;
                            } else {
                                break;
                            }
                        }
                }

                if (weightJ > 0) {
                    std::pair<int, int> edge = {I, weightJ};
                    /* should not need a single here since J is shared */
                    graph[J].push_back(edge);
                    // #pragma omp single
                    add_score_edge(scoreGraph, weightJ, J, I);
                    // std::cout << "New Edge: " << J << " -> " << I << " weightJ: " << weightJ << std::endl;
                }

                for (int i=fragments[I].length() -1, j=0;
                    i >=0 && j < fragments[J].length();
                    --i, ++j) {
                        weightI = 0;
                        for(int c = i, d = 0; c < fragments[I].length() && d < j; ++c, ++d) {
                            if(fragments[I][c] == fragments [J][d]) {
                                ++weightI;
                            } else {
                                break;
                            }
                        }
                }

                if (weightI > 0) {
                    std::pair<int, int> edge = {J, weightI};
                    // #pragma omp single
                    graph[I].push_back(edge);
                    // #pragma omp single
                    add_score_edge(scoreGraph, weightI, I, J);
                    // std::cout << "New Edge: " << I << " -> " << J << " weightI: " << weightI << std::endl;
                }

        }
    }
}

void print_adjacency_lists(graphLike & graph){
    std::cout << "Graph Adjacency Lists" << std::endl;

    for (auto it=graph.begin(); it!=graph.end(); ++it) {
        std::cout << it->first << ":\t";

        for (auto pair=it->second.begin(); pair!=it->second.end(); ++pair) {
            std::cout << pair->first << "(" << pair->second << "),\t";
        }
        std::cout << "degree: " << it->second.size();

        std::cout << std::endl;
    }
}

void merge_nodes(const int score, const int S, const int P, graphLike & graph, std::string & result) {
    //todo: somehow smash nodes into eachother? do i need reverse edges to do this effectivly?

}

void append_prefix_string(std::string & result, std:: string & P, int overlap) {
    std::cout << "Overlap: " << overlap << std::endl;
    std::cout << "Concatinating " << P << " with " << result << std::endl;
    result = P.substr(0, overlap - 1) + result;
    std::cout << "Result: " << result << std::endl;
}

void append_suffix_string(std::string & result, std::string & S, int overlap) {
    std::cout << "Overlap: " << overlap << std::endl;
    std::cout << "Concatinating " << result << " with " << S << std::endl;
    result += S.substr(overlap, S.length());
    std::cout << "Result: " << result << std::endl;
}


/* Greedily collapse the graph to a single node */
std::string colapse_graph(graphLike & graph, graphLike & scoreGraph, std::vector<std::string> & fragments) {
    std::string result;
    std::map<int, int> visited;
    for(int i=0; i < graph.size(); ++i) {
        visited[i] = 0;
    }

    int overlap = scoreGraph.rbegin()->first;
    std::pair<int, int> start_edge = scoreGraph.rbegin()->second[0];
    std::cout << "start node: " << fragments[start_edge.first] << std::endl;
    std::cout << "next node: " << fragments[start_edge.second] << std::endl;
    result = fragments[start_edge.first];
    visited[start_edge.first] = 1;

    append_suffix_string(result, fragments[start_edge.second], overlap);
    for (auto rit=scoreGraph.rbegin(); rit!=scoreGraph.rend(); ++rit) {
        std::cout << "Score: " << rit->first << " count: " << rit->second.size() << std::endl;
        //TODO: ...
    }

    return result;
}

int main(int argc, char**argv) {

    if (argc == 2) {
        std::vector<std::string> fragments;
        std::vector<unsigned long> fragment_lengths;
        unsigned long*packed_fragment;
		unsigned long num_blocks_in_fragment;
        std::vector<unsigned long*> packed_fragments;
		std::vector<unsigned long> num_blocks_in_fragments;

        graphLike graph;
        graphLike scoresGraph; //std::map<score, std::vector <std::pair<to, from>>>

        std::string packed_result;
        std::string result;

        
        /* read fragments */
		std::string fragment;
		std::ifstream fragments_fin(argv[1]);
		while(std::getline(fragments_fin, fragment)) {
			fragments.push_back(fragment);
			fragment_lengths.push_back(fragment.size());
		}
		fragments_fin.close();

        //TODO: maybe? idk
        // sort(fragments.begin(), fragments.end());

        /* bit pack fragments */
        // for(std::string fragment: fragments) {
		// 	bit_pack(packed_fragment, num_blocks_in_fragment, fragment);
		// 	packed_fragments.push_back(packed_fragment);
		// 	num_blocks_in_fragments.push_back(num_blocks_in_fragment);
		// }

        /* Sanity Check */
        // for(unsigned long i=0; i < fragments.size(); ++i) {
        //     std::cout << "original: " << fragments[i] << std::endl;
            // print_bitpacked_string(packed_fragments[i], num_blocks_in_fragments[i]);
            
            // std::string unpacked;
            // unpacked = unpack(packed_fragments[i], num_blocks_in_fragments[i], fragment_lengths[i]);
            // std::cout << "unpacked: " << unpacked << std::endl << std::endl;
        // }

        /* Build Graph */
        build_graph(fragments, graph, scoresGraph);

        // print_adjacency_lists(graph);
        
        /* Assemble... */
        result = colapse_graph(graph, scoresGraph, fragments);

        // std::string a = "ACT";
        // std::string b = "CTG";
        // int overlaping = 2;
        // append_suffix_string(a,b, overlaping);

        // b = "GAC";
        // append_prefix_string(a,b,overlaping);


        /* Output */
        // result = "GATTACCAATTACCAGGA"; /* TODO: DOTHIS" */
        std::cout << result << std::endl;
        return 0;
    }
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
}
