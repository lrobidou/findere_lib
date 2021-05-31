#pragma once
#include <fstream>
#include <string>
#include <zstr.hpp>

#include "utils.hpp"
inline std::string extractContentFromFasta(std::string filename) {
    std::ifstream myfilegz(filename);
    zstr::istream myfile(myfilegz);
    std::string line;
    std::string content;

    while (std::getline(myfile, line)) {
        if ((line[0] != '>') && (line[0] != '#')) {
            content += line;
        }
    }
    return content;
}

inline std::string extractContentFromFastqGz(std::string filename) {
    std::ifstream myfilegz(filename);
    zstr::istream myfile(myfilegz);

    std::string line;
    std::string content;

    int i = 0;
    while (std::getline(myfile, line)) {
        if (i == 1) {
            dbg("accept: " + line);
            content += line;
        } else {
            dbg("reject: " + line);
        }
        i++;
        i = i % 4;
    }

    return content;
}

inline std::string extractContentFromText(std::string filename) {
    std::ifstream myfilegz(filename);
    zstr::istream myfile(myfilegz);

    std::string line;
    std::string content;

    while (std::getline(myfile, line)) {
        content += line;
    }

    return content;
}