/**
 * Zachary Falkner
 * CSCI558 - Computational Biology
 * Cluster MS spectra
 * 
 * 
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


std::vector<std::vector<std::vector<int>>> locality_sensitive_hasing(std::vector<std::map<float, float>> & binned_spectras) {
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
    // min_y = -1 * max_y;

    min_x = *all_bins.begin();
    max_x = *all_bins.rbegin();

    // std::cout << "min_x: " << min_x << std::endl;
    // std::cout << "max_x: " << max_x << std::endl;
    // std::cout << "min_y: " << min_y << std::endl;
    // std::cout << "max_y: " << max_y << std::endl;

    //generate random points as pairs
    std::vector<std::pair<float, float>> points;
    for (size_t i=0; i<HASHES; ++i) {
        float x = random_float_in_range(min_x, max_x);
        float y = random_float_in_range(min_y, max_y);
        points.push_back({x,y});
    }

    std::vector<std::vector<std::vector<int>>> all_spectra_results;
    //for each spectra
    for (auto spectra: binned_spectras) {
        //for each random 2d point (hash, ~10)
        std::vector<std::vector<int>> hash_results;
        #pragma omp parallel for
        for (size_t p=0; p<HASHES; ++p) {
            //for each bin dot intensity with random 2d point to get angle
            std::vector<int> result;
            for (auto bin: spectra) {
                //TODO: audit this, not getting any negative numbers? maybe it's just the data
                float theta;
                theta = (points[p].first * bin.first) + (points[p].second * bin.second); /* A dot B */
                theta = theta / ( /* divide by magnitude of A mult magnitude of B */
                    sqrtf( (points[p].first * points[p].first) + (points[p].second * points[p].second) )
                    * sqrtf( (bin.first * bin.first) + (bin.second * bin.second) )
                );
                theta = acos(theta); /* arc cosis */
                theta = theta * (180.00f / PI);

                if (std::isnan(theta) == 1) {
                    theta = 0.0f;//don't break things
                }
                #pragma omp critical
                std::cout << "theta: " << theta << std::endl;

                //angle > 0 hash =1, angle < 0 hash =0
                int r;
                if (theta >= 0.0) {
                    r = 1;
                } else {
                    r = 0;
                }
                result.push_back(r);
            }
            #pragma omp critical
            hash_results.push_back(result);
        }
        all_spectra_results.push_back(hash_results);
    }
    return all_spectra_results;
}


/**
 * this is a very bad and insane approach i took, it wasn't working out,
 * essentially generating arbitrarily large normalized hyperplanes and measuring how they disect the universe of spectras
 * i was getting results, but i don't think they were good or useful
**/

// std::vector<std::vector<float>> locality_sensitive_hasing_other(std::vector<std::map<float, float>> & binned_spectras, std::set<float> & all_bins) {
//     //prints a nicely aligned picture of occupied bins
//     std::cout << std::setprecision(2) << std::fixed;
//     for (auto bin: all_bins) {
//         std::cout << bin << " ";
//     }
//     std::cout << std::endl;
//     for (auto spectra: binned_spectras) {
//         for (auto bin: all_bins) {
//             if(spectra.count(bin) != 0) {
//                 std::cout << bin << " ";
//             } else {
//                 std::cout << "       ";
//             }
//         }
//         std::cout << std::endl;
//     }
//     //end utility print


//     std::vector<std::map<float, float>> hashes;
//     std::vector<std::vector<float>> results;

//     // generate hashes, normalized hyper planes
//     for (size_t i=0; i<HASHES; ++i) {
//         std::map<float, float> hash;
//         size_t dimensions = all_bins.size();
//         for (auto bin: all_bins) {
//             int r = rand()%2;
//             hash[bin] = (1.0f/dimensions) * r;
//             // hash[bin] = r;
//         }
//         hashes.push_back(hash);
//     }

//     // normalize spectras
//     std::vector<float> magnitudes;
//     for (auto spectra: binned_spectras) {
//         float magnitude = 0.0;
//         for (auto it=spectra.begin(); it!=spectra.end(); ++it) {
//             magnitude += (it->second * it->second);
//         }
//         magnitude = sqrt(magnitude);
//         magnitudes.push_back(magnitude);
//     }
//     for (size_t i=0; i<binned_spectras.size(); ++i) {
//         auto spectra = binned_spectras[i];
//         auto magnitude = magnitudes[i];
//         for (auto it=spectra.begin(); it!=spectra.end(); ++it) {
//             it->second = (1.0f/magnitude) * it->second;
//         }
//     }

//     // on every spectra, apply every hash
//     for (auto spectra: binned_spectras) {
//         std::vector<float> spectra_results;
//         std::cout << "#spectra#" << std::endl;
//         for (auto hash: hashes) {
//             float result = 0.0;
//             for (auto it=spectra.begin(); it!=spectra.end(); ++it) {
//                 result += (it->second * hash[it->first]);
//             }
//             result /= (hash.size() * spectra.size());
//             result = acos(result) * 180.00 / PI;
//             spectra_results.push_back(result);
//             std::cout << "result: " << result << std::endl;
//         }
//         results.push_back(spectra_results);
//     }
//     return results;
// }

std::map<int, std::vector<int>> cluster(std::vector<std::vector<std::vector<int>>> & results) {
    std::map<int, std::vector<int>> clusters;
    // for (size_t i=0; i<results.size(); ++i) {
    //     for (size_t j=i+1; j<results.size(); ++j) {
    //         for (size_t k=0;k<results[i].size(); ++k) {
    //             float diff = fabs(results[i][k] - results[j][k]);
    //             std::cout << "k = " << k << " dif = " << diff << std::endl;
    //         }
    //     }
    // }

    for (auto spectra_results: results) {
        std::vector<unsigned int> keys = std::vector<unsigned int>(spectra_results[0].size());
        for (auto hash_result: spectra_results) {
            for (size_t r=0; r< hash_result.size(); ++r) {
                keys[r] <<= 1;
                keys[r] &= hash_result[r];
            }
        }
    }


}

void print_result_clusters(std::map<int, std::vector<int>> & results) {
    for (std::pair<int, std::vector<int>> result: results) {
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
        std::map<int, std::vector<int>> clusters;

        read_mgf_file(binned_spectras, argv[1]);

        lsh_results = locality_sensitive_hasing(binned_spectras);
        std::cout << "lsh complete" << std::endl;

        clusters = cluster(lsh_results);
        std::cout << "clustering done" << std::endl;
        
        print_result_clusters(clusters);
        return 0;
    }
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
}