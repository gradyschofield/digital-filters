//
// Created by Grady Schofield on 10/9/21.
//

#ifndef CPP_UTIL_H
#define CPP_UTIL_H

#include<vector>

using namespace std;

namespace Util {
    template<typename T>
    vector<T> shiftUp(vector<T> x) {
        vector<T> ret(1 + x.size());
        for(int i = 0; i < x.size(); ++i) {
            ret[i+1] = x[i];
        }
        return ret;
    }

linspace
argmin
objective
}

#endif //CPP_UTIL_H
