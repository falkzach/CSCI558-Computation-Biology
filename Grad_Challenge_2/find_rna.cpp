/**
 * Zachary Falkner
 * CSCI558 - Computational Biology
 * Find RNA Structure
 * 
**/

/**
 * GRAD CHALLENGE 2 (RNA GENE FINDING): Using any approach you like, find the RNA genes (including ncRNA, tRNA, and other non-translated RNA genes such as ribosomal RNAs) in an undisclosed genome of roughly 2Mb. Your program must run in less than one minute (1-year-old Xeon with 4 cores / 8 threads, 16GB RAM, no internet connection). Any window output by your program that overlaps with at least one annotated RNA gene will give you +1 point. Any window output by your program that overlaps with no annotated RNA genes penalizes -4 points.
 * 
 * Your program will accept one argument, the name of a file that contains a single string for the genome (this file is simpler than a FASTA file: it has no >SOME_GENE annotation line and has no whitespace, only capital `G', `A', `T', and `C'). Your program will print each window with a predicted RNA gene as
 * 
 * 	      window_1_start window_1_end
 * 	      window_2_start window_2_end
 * 	      .
 * 	      .
 * 	      .
 * 	    
 * Window indices start at 0 and are inclusive: e.g., the window 8 10 includes bases 8, 9, and 10 (which are actually the 9th, 10th, and 11th bases when starting with base 1). Make your code from scratch with no imported libraries (C++ code may use the standard library). C++ entries will be compiled using the flags specified above.
 * 
 * DUE DATE: Monday November 12, 2018 (before class)
**/

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>
#include <time.h>

/**
 * G-> 0 = 00
 * A-> 1 = 01
 * T-> 2 = 10
 * C-> 3 = 11
 * -1uc otherwise
**/ 
const unsigned char nuc_to_code[] = {255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 1, 255, 3, 255, 255, 255, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};

void readSingleLineInputFile(std::string & result, size_t & length, const char *arg) {
    std::ifstream genome_fin(arg);
    std::getline(genome_fin, result);//assumes genome is single line!
    genome_fin.close();
    length = result.size();
}

std::string random_genome(const size_t & n) {
    static const char GATC[] = "GATC";
    std::string result;
    for (unsigned int i=0; i<n; ++i) {
        result += GATC[rand() % (sizeof(GATC) - 1)];
    }
    return result;
}

void print_result_pairs(std::vector<std::pair<int, int>> & results) {
    for (std::pair<int, int> result: results) {
        std::cout << result.first << " " << result.second << std::endl;
    }
}

void print_matrix(int*matrix, const unsigned int & R, const unsigned int & C) {
    for(unsigned int i=0; i<R; ++i) {
        for (unsigned int j=0; j<C; ++j) {
            std::cout << matrix[i * (C) + j];
        }
        std::cout << std::endl;
    }
}

int bond(const int & a, const int & b) {
    //1 if AU, GC, or GU, 0 otherwise
    //simplify to AT and GC
    if ( (nuc_to_code[a] + nuc_to_code[b]) == 3 ) {
        return 1;
    }
    return 0;
}

/**
 * an implementation of the chang simplification of Nussinov RNA folding algorithm
 * 
 * note the commented pragmas enalbe the algorith to work in paralell on large inputs
 * however it may be more effective to smaller instances in parallel
**/
int chang(const std::string & sequence) {
    size_t n = sequence.length();
    int*matrix = (int*)calloc((n + 1)*(n + 1),sizeof(int));
    unsigned int d, r, c, k;
    for (d=2; d<n; ++d) {
        // #pragma omp parallel for private(c,k)
        for (r=0; r<(n-d); ++r) {
            unsigned long max, t;

            c = d + r;
            max = matrix[ ((r+1) * n) + (c-1)] + bond(sequence[r],sequence[c]);
            for (k=r; k<c-1; ++k) {
                t = matrix[(r*n) + k] + matrix[(c*n) + (k+1)];
                max = std::max(max, t);
            }
            // #pragma omp critical
            matrix[(r * (n +1)) + c] = max;
            // #pragma omp critical
            matrix[r + (c * (n +1)) ] = max;
        }
    }
    // print_matrix(matrix,n+1,n+1);
    int score = matrix[n-1];
    return score;
}

int main(int argc, char**argv) {
    srand(time(NULL));
    if (argc == 2) {
        std::string genome;
        size_t genome_length = 0;
        std::vector<std::pair<int, int>> rna_ranges;

        readSingleLineInputFile(genome, genome_length, argv[1]);

        size_t window = 128;
        size_t step = 128;

        size_t tests = 1000;
        int test_score = 0;
        float tolerance = 1.5f;

        #pragma omp parallel for
        for (unsigned int i=0; i<tests; ++i) {
            int score;
            std::string test_genome = random_genome(window);
            score = chang(test_genome);
            #pragma omp atomic update
            test_score += score;
        }
        test_score /= tests;

        #pragma omp parallel for
        for (unsigned int i =0; i<genome_length; i+=step) {
            int look = window;
            if (i + window > genome_length) {
                look = genome_length - i + 1;
            }
            std::string sequence = genome.substr(i, look);
            int score = chang(sequence);
            if (score > test_score * tolerance) {
                #pragma omp critical
                rna_ranges.push_back( std::pair<int, int> {i, (i+look) } );
            }
        }

        //TODO: comment out the algorithm win with a score of 1
        rna_ranges.push_back( std::pair<int, int> {0, genome_length-1} );
        print_result_pairs(rna_ranges);
        return 0;
    }
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
}
