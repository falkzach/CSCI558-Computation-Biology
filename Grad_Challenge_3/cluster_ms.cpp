/**
 * Zachary Falkner
 * CSCI558 - Computational Biology
 * Cluster MS spectra
 * 
 * This implementation attempts to use Locality Specific Hashing to Cluster MS Spectra.
 * I feel there are some flaws in my implementation leading to somewhat inconsistent clustering.
**/

/**
 *GRAD CHALLENGE 3 (CLUSTERING MS SPECTRA): Read in a .mgf file (example here) containing many thousand MS/MS spectra. Look only at the fragmentation peaks (not any precursor information). Discretize into m/z bins (by rounding) of 0.01 Th and cluster together the mass spectra. Your goals are (1) for each cluster has spectra that are not far from one another (using L2 distance) and (2) for spectra in different clusters are far from one another (using L2 distance).
 *
 *Your program will accept one argument, the name of the .mgf file. Your program will print the clusters by referring to the indices of the spectra in the .mgf file on the same line. For example, the output
 *
 *	      0 1 3 8
 *	      2 4 5
 *	      6 7
 *	      9
 *	    
 * has clusters {0,1,3,8}, {2,4,5}, {6,7}, and {9}. Distance between spectrum 0 and spectrum 1 will penalize this result, while distance between spectrum 0 and spectrum 2 will reward this result. Spectra unassigned to a cluster will be put into cluster 0. The winner will be selected by an F-test. Your code must run in <10 minutes on several thousand spectra (with # peaks roughly equal to one of the example spectra above). Deviate from the output format at your own risk.
 *
 * Make your code from scratch with no imported libraries (exception: C++ code may use the standard library and python may use numpy). C++ entries will be compiled using the flags specified above.
 *
 *DUE DATE: Monday December 3, 2018 (before class)
**/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include <math.h>
#include <omp.h>
#include <time.h>

#define PI 3.14159265

#define HASHES 10
#define BEGIN_IONS "BEGIN IONS"
#define END_IONS "END IONS"
constexpr char  READ_REGULAR_EXPRESSION[] = "^(\\d+\\.\\d+)\\s(\\d+\\.\\d+)$";

float random_float_in_range(float lo, float hi) {
    float r = lo + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(hi-lo)));
    return r;
}

/*
 * read in the mgf data binning it to 2 decimal places
 */
void read_mgf_file(std::vector<std::map<float, float>> & binned_spectras, const std::string & path) {
    std::string line;
    std::ifstream mgf_fin(path);
    bool inspectra = false;
    std::regex readRegex = std::regex(READ_REGULAR_EXPRESSION);
    std::smatch results;
    std::map<float, float> spectra;
    while(std::getline(mgf_fin, line)) {
        if(inspectra == true) {
            if (std::regex_search(line, results, readRegex)) {
                float massToCharge = roundf(std::atof(results[1].str().c_str()) * 100) / 100;  //rounds to 2 digits for binning
                float abundance = std::atof(results[2].str().c_str());
                if (spectra.count(massToCharge) == 0) {
                    spectra[massToCharge] = 0;
                }
                spectra[massToCharge] += abundance;
            } else if (line.compare(END_IONS) == 0) {
                binned_spectras.push_back(spectra);
                inspectra = false;
            }
        } else {
            if (line.compare(BEGIN_IONS) == 0) {
                spectra.clear();
                inspectra = true;
            }
        }
    }
}


std::vector<std::vector<std::vector<int>>> locality_sensitive_hasing(const std::vector<std::map<float, float>> & binned_spectras) {
    std::set<float> all_bins;
    float min_x = 0.0;
    float max_x = 0.0;
    float min_y = 0.0;
    float max_y = 0.0;

    for (auto spectra: binned_spectras) {
        for (auto bin: spectra) {
            if (all_bins.count(bin.first) == 0) {
                all_bins.insert(bin.first);
            }
            max_y = fmax(bin.second, max_y);
        }
    }

    min_x = *all_bins.begin();
    max_x = *all_bins.rbegin();

    //generate random points as pairs
    std::vector<std::pair<float, float>> points;
    for (size_t i=0; i<HASHES; ++i) {
        float x = random_float_in_range(min_x, max_x);
        float y = random_float_in_range(min_y, max_y);
        points.push_back({x,y});
    }

    std::vector<std::vector<std::vector<int>>> all_spectra_results = std::vector<std::vector<std::vector<int>>>(binned_spectras.size());
    //for each spectra
    #pragma omp parallel for
    for (size_t s=0; s<binned_spectras.size(); ++s) {
        auto spectra = binned_spectras[s];
        //for each random 2d point (hash, ~10)
        std::vector<std::vector<int>> hash_results;
        for (size_t p=0; p<HASHES; ++p) {
            //for each bin dot intensity with random 2d point to get angle
            std::vector<int> result;
            for (auto bin: spectra) {
                std::pair<float, float> x_axis = {1.0,0.0};
                float theta_1;
                theta_1 = (points[p].first * x_axis.first) + (points[p].second * x_axis.second); /* A dot B */
                theta_1 = theta_1 / ( /* divide by magnitude of A mult magnitude of B */
                    sqrtf( (points[p].first * points[p].first) + (points[p].second * points[p].second) )
                    * sqrtf( (x_axis.first * x_axis.first) + (x_axis.second * x_axis.second) )
                );
                theta_1 = acos(theta_1); /* arc cosis */
                theta_1 = theta_1 * (180.00f / PI);

                if (std::isnan(theta_1) == 1) {
                    theta_1 = 0.0f;//don't break things
                }

                float theta_2;
                theta_2 = (x_axis.first * bin.first) + (x_axis.second * bin.second); /* A dot B */
                theta_2 = theta_2 / ( /* divide by magnitude of A mult magnitude of B */
                    sqrtf( (x_axis.first * x_axis.first) + (x_axis.second * x_axis.second) )
                    * sqrtf( (bin.first * bin.first) + (bin.second * bin.second) )
                );
                theta_2 = acos(theta_2); /* arc cosis */
                theta_2 = theta_2 * (180.00f / PI);

                if (std::isnan(theta_2) == 1) {
                    theta_2 = 0.0f;//don't break things
                }

                int r;
                if (theta_2 >= theta_1) {
                    r = 1;
                } else {
                    r = 0;
                }
                result.push_back(r);
            }
            hash_results.push_back(result);
        }
        #pragma omp critical
        all_spectra_results[s] = hash_results;
    }
    return all_spectra_results;
}

std::map<int, std::set<int>> cluster(const std::vector<std::vector<std::vector<int>>> & results, const std::vector<std::map<float, float>> & binned_spectras) {
    std::map<int, std::set<int>> clusters;
    std::vector<int> cantor_hashes;

    /**
     * bit-pack the hash results of each bin together
     * TODO: this should be doable in the original scoring phase for later optimizations
    **/
    std::vector<std::vector<unsigned int>> spectra_keysets;
    for (auto spectra_results: results) {
        std::vector<unsigned int> keys = std::vector<unsigned int>(spectra_results[0].size(), 0);//ensure a 0 in each location
        for (auto hash_result: spectra_results) {
            for (size_t r=0; r< hash_result.size(); ++r) {
                keys[r] <<= 1;//shift in a 0
                keys[r] |= hash_result[r];//and in the result of the hash
            }
        }
        spectra_keysets.push_back(keys);
    }

    //cantor pairing, Robin Lockwood turned me onto this after much time spent fretting about "HOW DO I COMBINE ALL THESE HASHES"
    //TODO: don't think I'm using this correctly, clustering results are inconsistent!
    for (auto keyset: spectra_keysets) {
        unsigned long a = 0;
        unsigned long b = 0;
        size_t n = keyset.size();
        n = 2 << n;

        for (auto& n: keyset) {
            a += n;
            b += n;
        }
        b += 1;
        a /= n;
        b /= n;

        unsigned long h = a * b + 1;
        cantor_hashes.push_back(h);
    }

    for (size_t i=0; i<cantor_hashes.size(); ++i) {
        clusters[cantor_hashes[i]].insert(i);
    }
    return clusters;
}

void print_result_clusters(const std::map<int, std::set<int>> & results) {
    for (std::pair<int, std::set<int>> result: results) {
        for (int index: result.second) {
            std::cout << index << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char**argv) {
    srand(time(NULL));
    if (argc == 2) {
        std::vector<std::map<float, float>> binned_spectras;
        std::vector<std::vector<std::vector<int>>> lsh_results;
        std::map<int, std::set<int>> clusters;

        read_mgf_file(binned_spectras, argv[1]);
        lsh_results = locality_sensitive_hasing(binned_spectras);

        clusters = cluster(lsh_results, binned_spectras);
        
        print_result_clusters(clusters);
        return 0;
    }
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
}
