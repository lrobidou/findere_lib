#pragma once

#include <robin_hood.h>

#include <bf/all.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <zstr.hpp>

#include "../../src/reader.hpp"
#include "../../src/utils.hpp"

inline std::string make_canonical(std::string kmer) {
    std::string reverseComplement(kmer.rbegin(), kmer.rend());
    int l = kmer.length();
    for (int i = 0; i < l; i++) {
        switch (reverseComplement[i]) {
            case 'A':
                reverseComplement[i] = 'T';
                break;
            case 'C':
                reverseComplement[i] = 'G';
                break;
            case 'T':
                reverseComplement[i] = 'A';
                break;
            case 'G':
                reverseComplement[i] = 'C';
                break;
            case 'N':
                reverseComplement[i] = 'N';
                break;
            default:
                std::cerr << reverseComplement[i] << "encountered in a kmer. this is not alloex. Exiting." << std::endl;
                exit(1);
                break;
        }
        if (reverseComplement[i] > kmer[i]) {
            return kmer;
        }
    }
    return reverseComplement;
}

namespace truth {
inline robin_hood::unordered_set<std::string> indexFastas(std::vector<std::string> filenames, int k, bool canonical = false) {
    robin_hood::unordered_set<std::string> output;
    std::string line;
    for (auto const& filename : filenames) {
        std::ifstream myfile(filename);
        if (myfile.is_open()) {
            while (std::getline(myfile, line)) {
                if ((line[0] != '>') && (line[0] != '#')) {
                    unsigned long long start = 0;
                    unsigned long long l = line.length();
                    if (canonical) {
                        while ((start + k) <= l) {
                            output.insert(make_canonical(line.substr(start, k)));
                            start++;
                        }
                    } else {
                        while ((start + k) <= l) {
                            output.insert(line.substr(start, k));
                            start++;
                        }
                    }
                }
            }
            myfile.close();
        } else {
            std::cerr << "The file " << filename << " does not exist." << std::endl;
            exit(1);
        }
    }
    return output;
}

inline robin_hood::unordered_set<std::string> indexFastqGz(std::vector<std::string> filenames, int k, bool canonical = false) {
    robin_hood::unordered_set<std::string> output;
    std::string line;
    for (auto const& filename : filenames) {
        dbg("truth::indexFastqGz before declaring std::ifstream");
        std::ifstream myfilegz(filename);
        dbg("truth::indexFastqGz before declaring zstr::istream");
        zstr::istream myfile(myfilegz);
        dbg("truth::indexFastqGz starting iterating");
        int i = 0;
        while (std::getline(myfile, line)) {
            dbg("truth::indexFastqGz one new line");
            if (i == 1) {
                dbg("truth::indexFastqGz the line is relevant");
                unsigned long long start = 0;
                unsigned long long l = line.length();
                if (canonical) {
                    std::cout << "canonical" << std::endl;
                    while ((start + k) <= l) {
                        output.insert(make_canonical(line.substr(start, k)));
                        start++;
                    }
                } else {
                    while ((start + k) <= l) {
                        output.insert(line.substr(start, k));
                        start++;
                    }
                }

                dbg("truth::indexFastqGz the line is indexed");
            }
            i++;
            i = i % 4;
        }
    }
    return output;
}

inline robin_hood::unordered_set<std::string> indexText(std::vector<std::string> filenames, int k, bool canonical = false) {
    robin_hood::unordered_set<std::string> output;
    for (auto const& filename : filenames) {
        std::ifstream myfile(filename);
        std::string content = extractContentFromText(filename);
        unsigned long long start = 0;
        unsigned long long l = content.length();
        if (canonical) {
            while ((start + k) <= l) {
                output.insert(make_canonical(content.substr(start, k)));
                start++;
            }
        } else {
            while ((start + k) <= l) {
                output.insert(content.substr(start, k));
                start++;
            }
        }
    }

    return output;
}

}  // namespace truth

namespace findere_internal {
inline std::tuple<bf::basic_bloom_filter*, unsigned long long> indexFastasGivenBits(const std::vector<std::string>& filenames, unsigned long long bits, const unsigned numHashes, const unsigned int& k, const double& epsilon_percent, bool canonical = false) {
    bf::basic_bloom_filter* filter = new bf::basic_bloom_filter(bf::make_hasher(numHashes), bits);
    std::string line;
    for (auto const& filename : filenames) {
        std::ifstream myfile(filename);
        if (myfile.is_open()) {
            while (std::getline(myfile, line)) {
                if ((line[0] != '>') && (line[0] != '#')) {
                    unsigned long long start = 0;
                    unsigned long long l = line.length();
                    if (canonical) {
                        while ((start + k) <= l) {
                            filter->add(make_canonical(line.substr(start, k)));
                            start++;
                        }
                    } else {
                        while ((start + k) <= l) {
                            filter->add(line.substr(start, k));
                            start++;
                        }
                    }
                }
            }
            myfile.close();
        } else {
            std::cerr << "The file " << filename << " does not exist." << std::endl;
            exit(1);
        }
    }
    return {filter, bits};
}

inline std::tuple<bf::basic_bloom_filter*, unsigned long long> indexFastasGivenTruth(const std::vector<std::string>& filenames, const robin_hood::unordered_set<std::string>& truth, const unsigned numHashes, const unsigned int& k, const double& epsilon_percent, bool canonical = false) {
    // number of *unique* elements to add in that filter
    const unsigned long long n = truth.size();
    // size (in bit) required for that filter
    unsigned long long m = -(n / log(1 - ((double)epsilon_percent / (double)100)));
    // oops, maybe m is not a multiple of 8
    // this is required by most implementation of Bloom filters
    // let's fix that
    m = m + 8 - (m % 8);

    // now that we have the size, let's index those files
    return indexFastasGivenBits(filenames, m, numHashes, k, epsilon_percent, canonical);
}

inline std::tuple<bf::basic_bloom_filter*, unsigned long long> indexFastqGZGivenBits(const std::vector<std::string>& filenames, const unsigned long long bits, const unsigned numHashes, const unsigned int& k, bool canonical) {
    bf::basic_bloom_filter* filter = new bf::basic_bloom_filter(bf::make_hasher(numHashes), bits);

    std::string line;
    for (auto const& filename : filenames) {
        std::ifstream myfilegz(filename);
        zstr::istream myfile(myfilegz);
        int i = 0;
        while (std::getline(myfile, line)) {
            if (i == 1) {
                unsigned long long start = 0;
                unsigned long long l = line.length();
                if (canonical) {
                    while ((start + k) <= l) {
                        filter->add(make_canonical(line.substr(start, k)));
                        start++;
                    }
                } else {
                    while ((start + k) <= l) {
                        filter->add(line.substr(start, k));
                        start++;
                    }
                }
            }
            i++;
            i = i % 4;
        }
    }
    return {filter, bits};
}

inline std::tuple<bf::basic_bloom_filter*, unsigned long long> indexFastqGZGivenTruth(const std::vector<std::string>& filenames, const robin_hood::unordered_set<std::string>& truth, const unsigned numHashes, const unsigned int& k, const double& epsilon_percent, bool canonical) {
    // number of *unique* elements to add in that filter
    const unsigned long long n = truth.size();
    // size (in bit) required for that filter
    unsigned long long m = -(n / log(1 - ((double)epsilon_percent / (double)100)));
    // oops, maybe m is not a multiple of 8
    // this is required by most implementation of Bloom filters
    // let's fix that
    m = m + 8 - (m % 8);
    return indexFastqGZGivenBits(filenames, m, numHashes, k, canonical);
}

inline std::tuple<bf::basic_bloom_filter*, unsigned long long> indexTextGivenBits(const std::vector<std::string>& filenames, unsigned long long bits, const unsigned numHashes, const unsigned int& k, const double& epsilon_percent, bool canonical = false) {
    bf::basic_bloom_filter* filter = new bf::basic_bloom_filter(bf::make_hasher(numHashes), bits);

    for (auto const& filename : filenames) {
        // std::ifstream myfilegz(filename);
        // zstr::istream myfile(myfilegz);
        std::string content = extractContentFromText(filename);

        unsigned long long start = 0;
        unsigned long long l = content.length();
        if (canonical) {
            while ((start + k) <= l) {
                filter->add(make_canonical(content.substr(start, k)));
                start++;
            }
        } else {
            while ((start + k) <= l) {
                filter->add(content.substr(start, k));
                start++;
            }
        }
    }

    return {filter, bits};
}

inline std::tuple<bf::basic_bloom_filter*, unsigned long long> indexTextGivenTruth(const std::vector<std::string>& filenames, const robin_hood::unordered_set<std::string>& truth, const unsigned numHashes, const unsigned int& k, const double& epsilon_percent, bool canonical = false) {
    // number of *unique* elements to add in that filter
    const unsigned long long n = truth.size();
    // size (in bit) required for that filter
    unsigned long long m = -(n / log(1 - ((double)epsilon_percent / (double)100)));
    // oops, maybe m is not a multiple of 8
    // this is required by most implementation of Bloom filters
    // let's fix that
    m = m + 8 - (m % 8);
    return indexTextGivenBits(filenames, m, numHashes, k, epsilon_percent, canonical);
}
}  // namespace findere_internal

inline std::tuple<robin_hood::unordered_set<std::string>, bf::basic_bloom_filter*, int, unsigned long long> indexFastas(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, const double& epsilon_percent, bool canonical = false) {
    // create ground truth
    robin_hood::unordered_set<std::string> truth = truth::indexFastas(filenames, k, canonical);
    auto t0 = std::chrono::high_resolution_clock::now();
    const auto& [filter, sizeOfFilter] = findere_internal::indexFastasGivenTruth(filenames, truth, numHashes, k, epsilon_percent, canonical);
    auto t1 = std::chrono::high_resolution_clock::now();
    return {truth, filter, std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count(), sizeOfFilter};
}

inline std::tuple<robin_hood::unordered_set<std::string>, bf::basic_bloom_filter*, int, unsigned long long> indexFastqGz(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, const double& epsilon_percent, bool canonical = false) {
    // create ground truth
    robin_hood::unordered_set<std::string> truth = truth::indexFastqGz(filenames, k, canonical);
    auto t0 = std::chrono::high_resolution_clock::now();
    const auto& [filter, sizeOfFilter] = findere_internal::indexFastqGZGivenTruth(filenames, truth, numHashes, k, epsilon_percent, canonical);
    auto t1 = std::chrono::high_resolution_clock::now();
    return {truth, filter, std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count(), sizeOfFilter};
}

inline std::tuple<robin_hood::unordered_set<std::string>, bf::basic_bloom_filter*, int, unsigned long long> indexText(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, const double& epsilon_percent, bool canonical = false) {
    // create ground truth
    robin_hood::unordered_set<std::string> truth = truth::indexText(filenames, k, canonical);
    auto t0 = std::chrono::high_resolution_clock::now();
    const auto& [filter, sizeOfFilter] = findere_internal::indexTextGivenTruth(filenames, truth, numHashes, k, epsilon_percent, canonical);
    auto t1 = std::chrono::high_resolution_clock::now();
    return {truth, filter, std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count(), sizeOfFilter};
}

// passing bits as parameter

inline std::tuple<bf::basic_bloom_filter*, int, unsigned long long> indexFastasGivenBits(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, unsigned long long bits, bool canonical = false) {
    // create ground truth
    auto t0 = std::chrono::high_resolution_clock::now();
    const auto& [filter, sizeOfFilter] = findere_internal::indexFastasGivenBits(filenames, bits, numHashes, k, canonical);
    auto t1 = std::chrono::high_resolution_clock::now();
    return {filter, std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count(), sizeOfFilter};
}

inline std::tuple<bf::basic_bloom_filter*, int, unsigned long long> indexFastqGzGivenBits(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, unsigned long long bits, bool canonical = false) {
    // create ground truth
    auto t0 = std::chrono::high_resolution_clock::now();
    const auto& [filter, sizeOfFilter] = findere_internal::indexFastqGZGivenBits(filenames, bits, numHashes, k, canonical);
    auto t1 = std::chrono::high_resolution_clock::now();
    return {filter, std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count(), sizeOfFilter};
}

inline std::tuple<bf::basic_bloom_filter*, int, unsigned long long> indexTextGivenBits(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, unsigned long long bits, bool canonical = false) {
    // create ground truth
    auto t0 = std::chrono::high_resolution_clock::now();
    const auto& [filter, sizeOfFilter] = findere_internal::indexFastasGivenBits(filenames, bits, numHashes, k, canonical);
    auto t1 = std::chrono::high_resolution_clock::now();
    return {filter, std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count(), sizeOfFilter};
}

namespace findere {
inline std::tuple<robin_hood::unordered_set<std::string>, bf::basic_bloom_filter*, int, unsigned long long> indexFastas(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, const double& epsilon_percent, const unsigned& nbNeighboursMin, bool canonical = false) {
    // indexing for QTF is esay: just index as usual, but with k = k - nbNeighboursMin
    return ::indexFastas(filenames, numHashes, k - nbNeighboursMin, epsilon_percent, canonical);
}

inline std::tuple<robin_hood::unordered_set<std::string>, bf::basic_bloom_filter*, int, unsigned long long> indexFastqGz(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, const double& epsilon_percent, const unsigned& nbNeighboursMin, bool canonical = false) {
    // indexing for QTF is esay: just index as usual, but with k = k - nbNeighboursMin
    return ::indexFastqGz(filenames, numHashes, k - nbNeighboursMin, epsilon_percent, canonical);
}

inline std::tuple<robin_hood::unordered_set<std::string>, bf::basic_bloom_filter*, int, unsigned long long> indexText(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, const double& epsilon_percent, const unsigned& nbNeighboursMin, bool canonical = false) {
    // indexing for QTF is esay: just index as usual, but with k = k - nbNeighboursMin
    return ::indexText(filenames, numHashes, k - nbNeighboursMin, epsilon_percent, canonical);
}

// passing bits as parameter

inline std::tuple<bf::basic_bloom_filter*, int, unsigned long long> indexFastasGivenBits(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, const unsigned& nbNeighboursMin, unsigned long long bits, bool canonical = false) {
    // indexing for QTF is esay: just index as usual, but with k = k - nbNeighboursMin
    return ::indexFastasGivenBits(filenames, numHashes, k - nbNeighboursMin, bits, canonical);
}

inline std::tuple<bf::basic_bloom_filter*, int, unsigned long long> indexFastqGzGivenBits(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, const unsigned& nbNeighboursMin, unsigned long long bits, bool canonical = false) {
    // indexing for QTF is esay: just index as usual, but with k = k - nbNeighboursMin
    return ::indexFastqGzGivenBits(filenames, numHashes, k - nbNeighboursMin, bits, canonical);
}

inline std::tuple<bf::basic_bloom_filter*, int, unsigned long long> indexTextGivenBits(const std::vector<std::string>& filenames, const unsigned int& numHashes, const unsigned int& k, const unsigned& nbNeighboursMin, unsigned long long bits, bool canonical = false) {
    // indexing for QTF is esay: just index as usual, but with k = k - nbNeighboursMin
    return ::indexTextGivenBits(filenames, numHashes, k - nbNeighboursMin, bits, canonical);
}
}  // namespace findere