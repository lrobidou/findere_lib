#pragma once
#include <robin_hood.h>

#include <bf/all.hpp>

#include "../../src/utils.hpp"
#include "customAMQ.hpp"

inline bool oneQuery(bf::bloom_filter* filter, const std::string s) {
    return filter->lookup(s);
}

inline bool oneQuery(const robin_hood::unordered_set<std::string>& hashSet, const std::string& s) {
    return hashSet.contains(s);
}

inline bool oneQuery(const customAMQ& amq, const std::string& s) {
    return amq.contains(s);
}
template <typename T>
inline unsigned long long getNextPositiveKmerPositionInTheQuery(T filterOrTruth, const std::string& s, unsigned int k, const unsigned long long& nbNeighboursMin, unsigned long long j, unsigned long long& nbQuery) {
    // std::cout << "getNextPositiveKmerPositionInTheQuery: " << k << " " << j << std::endl;
    unsigned long long numberOfJumps = 0;
    unsigned long long size = s.size();
    do {
        // std::cout << j + numberOfJumps * nbNeighboursMin << " is false." << std::endl;
        numberOfJumps += 1;
        nbQuery++;
    } while ((j + numberOfJumps * nbNeighboursMin < size - k + 1) && (oneQuery(filterOrTruth, s.substr(j + numberOfJumps * nbNeighboursMin, k))) == false);
    numberOfJumps -= 1;
    // std::cout << "returning numberOfJumps = numberOfJumps" << numberOfJumps << std::endl;
    // if (!(j + nbNeighboursMin < size - k + 1)) {
    //     // no positive kmer next
    //     std::cout << "j + nbNeighboursMin if after the end of the query" << std::endl;
    //     return j;
    // }
    // std::cout << "j=" << j << " si a positive" << std::endl;

    // j is now the index of a posiive kmer in the query. Let's backtrack, one kmer by one kmer, until we fond a negative kmer.

    //TODO eviter dernier query
    // while (oneQuery(filterOrTruth, s.substr(j - 1, k))) {
    //     std::cout << j - 1 << " is a positive" << std::endl;
    //     j--;
    //     nbQuery++;
    // }

    // nbQuery++;
    return numberOfJumps;
}

namespace truth {
inline std::vector<bool> queryTruth(const robin_hood::unordered_set<std::string>& truth, const std::string& seq, int k) {
    unsigned long long start = 0;
    unsigned long long l = seq.length();
    std::vector<bool> response(l - k + 1);

    while ((start + k) <= l) {
        response[start] = truth.contains(seq.substr(start, k));
        start++;
    }
    return response;
}
}  // namespace truth

namespace noQTF {
inline std::vector<bool> query(bf::bloom_filter* filter, const std::string& s, const unsigned int& k) {
    checknonNull(filter, "Nullptr passed to query function.");
    unsigned long long size = s.size();
    std::vector<bool> response(size - k + 1);

    for (unsigned long long i = 0; i < size - k + 1; i++) {
        response[i] = filter->lookup(s.substr(i, k));
    }

    return response;
}
}  // namespace noQTF

//TODO move to internal
template <typename T>
inline std::vector<bool> queryFilterOrTruth(const T& filterOrTruth, const std::string& s, const unsigned int& k, const unsigned int& nbNeighboursMin, unsigned long long& nbStretch, bool skip = false) {
    const unsigned int smallK = k - nbNeighboursMin;
    unsigned long long size = s.size();
    std::vector<bool> response(size - smallK + 1 - nbNeighboursMin);
    unsigned long long i = 0;              // index of the response vector
    unsigned long long stretchLength = 0;  // number of consecutive positives kmers
    unsigned long long j = 0;              // index of the query vector

    while (j < size - k + 1) {
        if (oneQuery(filterOrTruth, s.substr(j, smallK))) {
            stretchLength++;
            j++;
        } else {
            if (stretchLength != 0) {
                if (stretchLength > nbNeighboursMin) {
                    nbStretch++;
                    for (unsigned long long t = 0; t < stretchLength - nbNeighboursMin; t++) {
                        response[i] = true;
                        i++;
                    }
                    for (unsigned long long t = 0; t < nbNeighboursMin; t++) {
                        response[i] = false;
                        i++;
                    }
                } else {
                    for (unsigned long long t = 0; t < stretchLength; t++) {
                        response[i] = false;
                        i++;
                    }
                }
                stretchLength = 0;
            }

            // skip queries between current position and the next positive kmer
            if (skip && (nbNeighboursMin > 0)) {
                unsigned long long dontCare = 0;
                unsigned long long numberOfJumps = getNextPositiveKmerPositionInTheQuery(filterOrTruth, s, smallK, nbNeighboursMin, j, dontCare);
                // std::cout << "numberOfJumps = " << numberOfJumps << std::endl;
                for (unsigned long long temp = 0; temp < nbNeighboursMin * numberOfJumps; temp++) {
                    response[i] = false;
                    i++;
                    j++;
                }
            }

            response[i] = 0;
            i++;
            j++;
        }
    }
    if (stretchLength != 0) {
        if (stretchLength > nbNeighboursMin) {
            nbStretch++;
            for (unsigned long long k = 0; k < stretchLength; k++) {
                response[i] = true;
                i++;
            }
        } else {
            for (unsigned long long k = 0; k < stretchLength; k++) {
                response[i] = false;
                i++;
            }
        }
        stretchLength = 0;
    }
    return response;
}

namespace findere {

std::vector<bool> inline query(bf::basic_bloom_filter* filter, const std::string& s, const unsigned int& k, const unsigned int& nbNeighboursMin, unsigned long long& nbStretch, bool skip = false) {
    return queryFilterOrTruth(filter, s, k, nbNeighboursMin, nbStretch, skip);
}

std::vector<bool> inline query(const robin_hood::unordered_set<std::string>& truth, const std::string& s, const unsigned int& k, const unsigned int& nbNeighboursMin, unsigned long long& nbStretch, bool skip = false) {
    return queryFilterOrTruth(truth, s, k, nbNeighboursMin, nbStretch, skip);
}

std::vector<bool> inline query(bf::basic_bloom_filter* filter, const std::string& s, const unsigned int& k, const unsigned long long& nbNeighboursMin, bool skip = false) {
    unsigned long long dontCare = 0;
    return queryFilterOrTruth(filter, s, k, nbNeighboursMin, dontCare, skip);
}

std::vector<bool> inline query(const robin_hood::unordered_set<std::string>& truth, const std::string& s, const unsigned int& k, const unsigned int& nbNeighboursMin, bool skip = false) {
    unsigned long long dontCare = 0;
    return queryFilterOrTruth(truth, s, k, nbNeighboursMin, dontCare, skip);
}

std::vector<bool> inline query(const customAMQ& amq, const std::string& s, const unsigned int& k, const unsigned int& nbNeighboursMin, unsigned long long& nbStretch, bool skip = false) {
    return queryFilterOrTruth(amq, s, k, nbNeighboursMin, nbStretch, skip);
}

std::vector<bool> inline query(const customAMQ& amq, const std::string& s, const unsigned int& k, const unsigned long long& nbNeighboursMin, bool skip = false) {
    unsigned long long dontCare = 0;
    return queryFilterOrTruth(amq, s, k, nbNeighboursMin, dontCare, skip);
}

}  // namespace findere