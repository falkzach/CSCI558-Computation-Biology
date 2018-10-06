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
        
        
        
        
        std::vector<unsigned long> fragment_lengths;
        unsigned long*packed_fragment;
		unsigned long num_blocks_in_fragment;
        std::vector<unsigned long*> packed_fragments;
		std::vector<unsigned long> num_blocks_in_fragments;        


        fragment_lengths.push_back(fragment.size());
        
        
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