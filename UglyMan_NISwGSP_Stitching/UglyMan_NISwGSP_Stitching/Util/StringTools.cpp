//
//  StringTools.cpp
//  UglyMan_NISwGSP_Stitching
//
//  Created by RiverHe on 17/10/16.
//  Copyright © 2017年 nothinglo. All rights reserved.
//

#include "StringTools.h"

vector<string> split(string str, string pattern) {
    vector<string> ret;
    if (pattern.empty()) {
        return ret;
    }
    size_t start = 0, index = str.find_first_of(pattern, 0);
    while (index != str.npos) {
        if (start != index) {
            ret.push_back(str.substr(start, index - start));
        }
        start = index + 1;
        index = str.find_first_of(pattern, start);
    }
    if (!str.substr(start).empty()) {
        ret.push_back(str.substr(start));
    }
    return ret;
}