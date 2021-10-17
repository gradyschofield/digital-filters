//
// Created by Grady Schofield on 10/16/21.
//

#ifndef CPP_BSPLINE_H
#define CPP_BSPLINE_H

#include<vector>

#include<Util.h>

using namespace std;

template<typename T>
class BSpline {
    int degree;
    vector<T> knots;
    vector<vector<T>> bases;

public:
    static vector<T> getDefaultKnots() {
        return vector<T>{0, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 48000};
    }

    int getNumBases() const {
        return bases.size();
    }

    vector<T> getBasis(int i) const {
        return bases[i];
    }

    static T evaluateSpline0(int i, vector<T> const & knots, T x) {
        return knots[i] <= x && x <= knots[i + 1] ? 1 : 0;
    }

    static T evaluateSpline(int i, int p, vector<T> const & knots, T x) {
        if(p == 0) {
            return BSpline::evaluateSpline0(i, knots, x);
        }
        T d1 = knots[i+p] - knots[i];
        T c1 = d1 == 0 ? 0 : (x - knots[i])/d1 * BSpline::evaluateSpline(i, p-1, knots, x);
        T d2 = knots[i+p+1] - knots[i+1];
        T c2 = d2 == 0 ? 0 : (knots[i+p+1] - x)/d2 * BSpline::evaluateSpline(i+1, p-1, knots, x);
        return c1 + c2;
    }

    BSpline(vector<T> const & grid, vector<T> const & knotsParam, int degree = 3)
        : degree(degree), knots(knotsParam)
    {
        for(int i = 0; i < degree; ++i) {
            knots.insert(begin(knots), knots[0]);
            knots.push_back(knots.back());
        }
        int numBases = knots.size()-degree-1;
        for(int i = 0; i < numBases; ++i) {
            auto basis = [&](T x){
                return BSpline::evaluateSpline(i, degree, knots, x);
            };
            vector<T> basisVector(grid.size());
            for(int j = 0; j < grid.size(); ++j) {
                basisVector[j] = grid[j] <= knots[i+degree+1] && grid[j] >= knots[i] ? basis(grid[j]) : 0;
            }
            bases.push_back(move(basisVector));
        }
    }

};

#endif //CPP_BSPLINE_H
