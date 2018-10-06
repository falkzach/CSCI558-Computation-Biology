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
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

/* graphLike is a map of node id's to vectors of pairs of node id and weight */
using graphLike = std::map<int, std::vector<std::pair<int, int>>>;

void prune_substring_fragments(std::vector<std::string> & fragments) {
    std::set<int> prune_fragments;
    for(int I=0; I < fragments.size(); ++I) {
        for (int J=I+1; J < fragments.size(); ++J) {
            if(fragments[I].length() <= fragments[J].length()) {
                if(fragments[J].find(fragments[I]) != std::string::npos) { // I is substring of J, erase I
                    prune_fragments.insert(I);
                }
            } else {
                if(fragments[I].find(fragments[J]) != std::string::npos) { // J is substring of I, erase J
                    prune_fragments.insert(J);
                }
            }
        }
    }
    for(std::set<int>::reverse_iterator rit=prune_fragments.rbegin(); rit != prune_fragments.rend(); rit++) {
        int id = *rit;
        fragments.erase(fragments.begin() + id);
    }

    // for(int I=0; I < fragments.size(); ++I) {
    //     std::cout<<fragments[I]<<std::endl;
    // }
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
void build_graph(const std::vector<std::string> fragments, graphLike & graph, graphLike & reverseGraph, graphLike & scoreGraph) {
    int min_overlap = 3;
    /* put all nodes in the graph with no edges */
    // #pragma omp parallel for
    for(unsigned long i=0; i < fragments.size(); ++i) {
        std::vector<std::pair<int, int>> adjacencies;
        std::vector<std::pair<int, int>> reverseAdjacencies;
        graph[i] = adjacencies;
        reverseGraph[i] = reverseAdjacencies;
    }

    unsigned long I, J;
    /* build out edges of graph */
    // #pragma omp parallel for private(J, weightI, weightJ, bestI, bestJ) shared(I)
    for(I=0; I < fragments.size(); ++I) {
        for (J=I+1; J < fragments.size(); ++J) {
            int weightI = 0;  //suffix of I prefix of J I -> J
            int weightJ = 0;  //suffix of J prefix of I J -> I
            int bestI = 0;
            int bestJ = 0;
            for (int i = 0, j=fragments[J].length() - 1; i < fragments[I].length() && j >= 0; ++i, --j) {
                    weightJ = 0;
                    for (int c = 0, d = j; c <= i && d < fragments[J].length(); ++c, ++d) {
                        if(fragments[I][c] == fragments [J][d]) {
                            ++weightJ;
                        } else {
                            weightJ = 0;
                            break;
                        }
                    }
                    bestJ = std::max(bestJ, weightJ);
            }

            for (int i=fragments[I].length() -1, j=0; i >=0 && j < fragments[J].length(); --i, ++j) {
                    weightI = 0;
                    for(int c = i, d = 0; c < fragments[I].length() && d <= j; ++c, ++d) {
                        if(fragments[I][c] == fragments [J][d]) {
                            ++weightI;
                        } else {
                            weightI = 0;
                            break;
                        }
                    }
                    bestI = std::max(bestI, weightI);
            }

            if (bestJ >= min_overlap) {
                std::pair<int, int> edge = {I, bestJ};
                std::pair<int, int> reverseEdge = {J, bestJ};
                graph[J].push_back(edge);
                reverseGraph[I].push_back(reverseEdge);
                add_score_edge(scoreGraph, bestJ, J, I);
            }

            if (bestI >= min_overlap) {
                std::pair<int, int> edge = {J, bestI};
                std::pair<int, int> reverseEdge = {I, bestI};
                graph[I].push_back(edge);
                reverseGraph[J].push_back(reverseEdge);
                add_score_edge(scoreGraph, bestI, I, J);
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
        // std::cout << "degree: " << it->second.size();

        std::cout << std::endl;
    }
}

void append_prefix_string(std::string & result, std:: string & P, int overlap) {
    result = P.substr(0, overlap) + result;
    std::cout << result << std::endl;
}

void append_suffix_string(std::string & result, std::string & S, int overlap) {
    result += S.substr(overlap + 1, S.length());
}

int all_set(std::map<int, int> visited) {
    int n_set = 0;
    n_set = std::accumulate(visited.begin(), visited.end(), 0, [] (int value, const std::map<int, int>::value_type& p) { return value + p.second; });
    return n_set == visited.size();
}

/* Greedily collapse the graph to a single node */
std::string colapse_graph(graphLike & graph, graphLike & scoreGraph, std::vector<std::string> & fragments) {
    std::string result;
    int result_head, result_tail;
    std::map<int, int> visited;
    for(int i=0; i < graph.size(); ++i) {
        visited[i] = 0;
    }

    // result = fragments[edge.first];
    for (auto score_rit=scoreGraph.rbegin(); score_rit!=scoreGraph.rend(); ++score_rit) {
        std::cout << "Score: " << score_rit->first << " count: " << score_rit->second.size() << std::endl;
        int overlap = score_rit->first;
        for(int i=0; i < score_rit->second.size(); ++i) {
            std::pair<int, int> edge = score_rit->second[i];
            visited[edge.first] = 1;
            if (all_set(visited) == 0) {
                edge  = score_rit->second[i];
                if(edge.second == result_head) {
                    std::cout << "a: " << fragments[edge.first] << std::endl;
                    std::cout << "r: " << result << std::endl;
                    append_prefix_string(result, fragments[edge.second], overlap);
                } else if(visited.count(edge.first) == 0) {
                    std::cout << "r: " << result << std::endl;
                    std::cout << "b: " << fragments[edge.second] << std::endl;
                    append_suffix_string(result, fragments[edge.second], overlap);
                }
            } else {
                break;
            }
        }
        if(all_set(visited)) { break; }
    }

    return result;
}

std::string walk_graph(graphLike & graph, graphLike & scoreGraph, std::vector<std::string> & fragments) {
    std::string result;
    std::map<int, int> path;
    std::map<int, int> overlaps;
    std::set<int> visited;
    int current_node;

    for (auto score_rit=scoreGraph.rbegin(); score_rit!=scoreGraph.rend(); ++score_rit) {
        int overlap = score_rit->first;
        for(int i=0; i < score_rit->second.size(); ++i) {
            std::pair<int, int> edge = score_rit->second[i];
            if(visited.count(edge.second) == 0  && path.count(edge.first) == 0) {
                visited.insert(edge.second);
                path[edge.first] = edge.second;
                overlaps[edge.first] = overlap;
            }
        }
    }

    // std::cout << "Path" << std::endl;
    // for(std::map<int,int>::iterator it=path.begin(); it!=path.end(); it++) {
    //     std::cout << it->first << " -> " << it->second << std::endl;
    // }

    
    for(int i=0; i< fragments.size(); ++i) {
        if(visited.count(i) == 0) {
            current_node = i;
            break;
        }
    }
    // std::cout << "Start Node: " << current_node << std::endl;
    result = fragments[current_node];
    std::cout << result << std::endl;
    int overlap; 
    int offset = result.length() - overlaps[current_node]; 
    while (path.count(current_node) > 0) {
        overlap = overlaps[current_node];
        for(int o=0; o<offset;++o) { std::cout << " "; }
        std::cout << fragments[path[current_node]] << std::endl;
        int append_node = current_node;
        current_node = path[current_node];
        offset += fragments[current_node].length() - overlaps[current_node];
        append_suffix_string(result, fragments[current_node], 0); //TODO: fix overlap lookups!!!
        // std::cout << result << std::endl;
    }

    return result;
}

int main(int argc, char**argv) {

    if (argc == 2) {
        std::vector<std::string> fragments;
        graphLike graph;
        graphLike reverseGraph;
        graphLike scoresGraph; //std::map<score, std::vector <std::pair<to, from>>>
        std::string result;

        /* read fragments */
        std::string fragment;
        std::ifstream fragments_fin(argv[1]);
        while(std::getline(fragments_fin, fragment)) {
            /* clean newlines and returns */
            fragment.erase(std::remove(fragment.begin(), fragment.end(), '\n'), fragment.end());
            fragment.erase(std::remove(fragment.begin(), fragment.end(), '\r'), fragment.end());
            fragments.push_back(fragment);
        }
        fragments_fin.close();

        /* prune pure substrings that add no new information */
        prune_substring_fragments(fragments);

        // std::cout<<"pruned"<<std::endl;
        // std::cout<<"n: " << fragments.size() << std::endl;

        /* Build Graph */
        build_graph(fragments, graph, reverseGraph, scoresGraph);

        // print_adjacency_lists(graph);
        // print_adjacency_lists(reverseGraph);
        
        /* Assemble... */
        // result = colapse_graph(graph, scoresGraph, fragments);
        result = walk_graph(graph, scoresGraph, fragments);

        /* Output */
        // result = "GATTACCAATTACCAGGA"; /* TODO: DOTHIS" */
        std::cout << result << std::endl;
        return 0;
    }
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
}
