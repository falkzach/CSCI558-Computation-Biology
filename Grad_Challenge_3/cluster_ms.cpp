/**
 * Zachary Falkner
 * CSCI558 - Computational Biology
 * Cluster MS Spectre
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
// #define READ_PATTERN = "^(\\d+\\.\\d+)\\s(\\d+\\.\\d+)$"


/*
 * http://fiehnlab.ucdavis.edu/projects/lipidblast/mgf-files
 */
void read_mgf_file(std::vector<std::map<float, float>> & binned_spectres, std::set<float> & all_bins, const std::string & path) {
    std::string line;
    std::ifstream mgf_fin(path);
    bool inSpectre = false;
    std::regex readRegex = std::regex("^(\\d+\\.\\d+)\\s(\\d+\\.\\d+)$");
    std::smatch results;
    std::map<float, float> spectre;
    while(std::getline(mgf_fin, line)) {
        if(inSpectre == true) {
            if (std::regex_search(line, results, readRegex)) {
                float massToCharge = roundf(std::atof(results[1].str().c_str()) * 100) / 100;  //rounds to 2 digits prior to binning
                float abundance = std::atof(results[2].str().c_str());
                if (spectre.count(massToCharge) == 0) {
                    spectre[massToCharge] = 0;
                }
                if (all_bins.count(massToCharge) == 0) {
                    all_bins.insert(massToCharge);
                }
                spectre[massToCharge] += abundance;
            } else if (line.compare(END_IONS) == 0) {
                binned_spectres.push_back(spectre);
                inSpectre = false;
            }
        } else {
            if (line.compare(BEGIN_IONS) == 0) {
                spectre.clear();
                inSpectre = true;
            }
        }
    }
}

void locality_sensitive_hasing(std::vector<std::map<float, float>> & binned_spectres, std::set<float> & all_bins) {
    //prints a nicely aligned picture of occupied bins
    std::cout << std::setprecision(2) << std::fixed;
    for (auto bin: all_bins) {
        std::cout << bin << " ";
    }
    std::cout << std::endl;
    for (auto spectre: binned_spectres) {
        for (auto bin: all_bins) {
            if(spectre.count(bin) != 0) {
                std::cout << bin << " ";
            } else {
                std::cout << "       ";
            }
        }
        std::cout << std::endl;
    }
    //end utility print


    std::vector<std::map<float, float>> hashes;
    std::vector<std::vector<float>> results;

    // generate hashes, normalized hyper planes
    for (size_t i=0; i<HASHES; ++i) {
        std::map<float, float> hash;
        size_t dimensions = all_bins.size();
        for (auto bin: all_bins) {
            int r = rand()%2;
            hash[bin] = (1.0f/dimensions) * r;
            // hash[bin] = r;
        }
        hashes.push_back(hash);
    }

    // normalize spectres
    std::vector<float> magnitudes;
    for (auto spectre: binned_spectres) {
        float magnitude = 0.0;
        for (auto it=spectre.begin(); it!=spectre.end(); ++it) {
            magnitude += (it->second * it->second);
        }
        magnitude = sqrt(magnitude);
        magnitudes.push_back(magnitude);
    }
    for (size_t i=0; i<binned_spectres.size(); ++i) {
        auto spectre = binned_spectres[i];
        auto magnitude = magnitudes[i];
        for (auto it=spectre.begin(); it!=spectre.end(); ++it) {
            it->second = (1.0f/magnitude) * it->second;
        }
    }

    // on every spectre, apply every hash
    for (auto spectre: binned_spectres) {
        std::vector<float> spectre_results;
        std::cout << "#spectre#" << std::endl;
        for (auto hash: hashes) {
            float result = 0.0;
            for (auto it=spectre.begin(); it!=spectre.end(); ++it) {
                result += (it->second * hash[it->first]);
            }
            result /= (hash.size() * spectre.size());
            result = acos(result) * 180.00 / PI;
            spectre_results.push_back(result);
            std::cout << "result: " << result << std::endl;
        }
        results.push_back(spectre_results);
    }



}

void cluster(std::vector<std::vector<float>> & results) {
    // for (size_t i=0; i<results.size(); ++i) {
    //     for (size_t j=i+1; j<results.size(); ++j) {
            
    //     }
    // }

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
        std::vector<std::map<float, float>> binned_spectres;
        std::set<float> all_bins;
        std::map<int, std::vector<int>> clusters;

        //TODO: implement....

        read_mgf_file(binned_spectres, all_bins, argv[1]);

        locality_sensitive_hasing(binned_spectres, all_bins);
        
        print_result_clusters(clusters);
        return 0;
    }
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
}