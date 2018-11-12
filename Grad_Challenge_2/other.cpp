
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
	std::cout << std::dec << std::endl;
}

std::string unpack_bit_packed_genome(const unsigned long * genome, const size_t & num_blocks, size_t & length) {
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

int get_packed_value_at_index(const unsigned long * genome, const size_t & index) {
    unsigned long result, offset;
    result = genome[(index / NUCLEOTITES_PER_BLOCK)];
    offset = index % NUCLEOTITES_PER_BLOCK;
    result <<= (2 * offset);
    result &= nibs_to_mask[1];
    result >>= (NUCLEOTITES_PER_BLOCK * 2) - 2;
    return result;
}


int bond_bit_packed(const int & a, const int & b) {
    //1 if AU, GC, or GU, 0 otherwise
    //simplify to AT and GC
    if ( (a + b) == 3 )
    {
        return 1;
    }
    return 0;
}

/**
 * https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-S8-S1
**/
int chang_bit_packed(const unsigned long * sequence, const size_t & num_blocks, const size_t & n, std::vector<std::pair<int,int>> & results) {
    int*matrix = (int*)calloc((n + 1)*(n + 1),sizeof(int));
    unsigned int d, r, c, k;
    for (d=2; d<n; ++d) {
        #pragma omp parallel for private(c,k)
        for (r=0; r<(n-d); ++r) {
            unsigned int a, b;
            long max, t;
            
            c = d + r;
            a = get_packed_value_at_index(sequence, r);
            b = get_packed_value_at_index(sequence, c);

            max = matrix[ ((r+1) * n) + (c-1)] + bond_bit_packed(a,b);
            for (k=r; k<c; ++k) {
                t = matrix[(r*n) + k] + matrix[(c*n) + (k+1)];
                max = std::max(max, t);
            }
            #pragma omp critical
            matrix[(r * (n +1)) + c] = max;
            #pragma omp critical
            matrix[(c * (n + 1)) + r] = max;
        }
    }
    // print_matrix(matrix, n + 1, n + 1);
    int score = matrix[n-1];
    return score;
}



// unsigned long * packed_genome;
// size_t num_blocks_in_genome;
// std::string unpacked_genome;



// bit_pack(packed_genome, num_blocks_in_genome, genome);

// if (SANITY_CHECK) {
//     std::cout<<"original: " << genome << std::endl;
//     print_bitpacked_string(packed_genome, num_blocks_in_genome);

//     unpacked_genome = unpack_bit_packed_genome(packed_genome, num_blocks_in_genome, genome_length);
//     std::cout << unpacked_genome << std::endl;
// }

//running chang on the bitpacked genome yields no benefit over running it on the string
// int score = chang_bit_packed(packed_genome, num_blocks_in_genome, genome_length, rna_ranges);