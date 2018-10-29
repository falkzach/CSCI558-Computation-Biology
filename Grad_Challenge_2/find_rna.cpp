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
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/**
 * TODO: disable debugging
 **/ 
#define SANITY_CHECK 1
#define DEBUG_ON 1

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

void readSingleLineInputFileFromArgument(std::string & result, size_t & length, const char *arg) {
        std::ifstream genome_fin(arg);
        std::getline(genome_fin, result);//assumes genome is single line!
        genome_fin.close();
        length = result.size();
}

void bit_pack(unsigned long * & result, size_t & num_blocks, const std::string & genome) {

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

			next_block |= nuc_to_code[ (int) genome[i] ];
		}
		result[block] = next_block;
	}

	// this loop is the inner body of the last itteration of the above loop
	unsigned long next_block = 0;
	for (unsigned char j=0; j<NUCLEOTITES_PER_BLOCK; ++j, ++i) {
		next_block <<= 2;

		if (i < genome.size()) {
			next_block |= nuc_to_code[ (int) genome[i] ];
		}
	}
	result[num_blocks-1] = next_block;
}

void print_bitpacked_string(const unsigned long* genome, const size_t & num_blocks) {
    std::cout << "packed: ";
	for (unsigned long block=0; block<num_blocks; ++block) {
		std::cout << std::hex << genome[block];
	}
	std::cout << std::endl;
}

std::string unpack_bit_packed_genome(const unsigned long * genome, const size_t & num_blocks, size_t length) {
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

void print_result_pairs(std::vector<std::pair<int, int>> & results) {
    for (std::pair<int, int> result: results) {
        std::cout << result.first << " " << result.second << std::endl;
    }

}

int main(int argc, char**argv) {
    
    if (argc == 2) {

        std::string genome;
        size_t genome_length = 0;
        unsigned long * packed_genome;
        size_t num_blocks_in_genome;
        std::string unpacked_genome;

        std::vector<std::pair<int, int>> rna_ranges;

        readSingleLineInputFileFromArgument(genome, genome_length, argv[1]);
        bit_pack(packed_genome, num_blocks_in_genome, genome);

        if (SANITY_CHECK) {
            std::cout<<"original: " << genome << std::endl;
            print_bitpacked_string(packed_genome, num_blocks_in_genome);

            unpacked_genome = unpack_bit_packed_genome(packed_genome, num_blocks_in_genome, genome_length);
            std::cout << unpacked_genome << std::endl;
        }

        // rna_ranges.push_back( std::pair<int, int> {1,2} );
        print_result_pairs(rna_ranges);

        return 0;
    }
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
}
