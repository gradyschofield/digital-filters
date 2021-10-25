//
// Created by Grady Schofield on 10/16/21.
//

#ifndef CPP_BSPLINE_H
#define CPP_BSPLINE_H

#include<vector>

#include<Util.h>
#include<DerivativeOptimizable.h>
#include<GeneticallyOptimizable.h>

using namespace std;

template<typename T>
class BSpline : public DerivativeOptimizable<BSpline, T>, public GeneticallyOptimizable<BSpline, T> {
    int degree;
    vector<T> knots;

public:
    static vector<T> getDefaultKnots() {
        return vector<T>{0, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 48000};
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

    static int getNumBases(int numKnots, int degree) {
        return numKnots-degree-1;
    }

    BSpline(vector<T> const & knotsParam, int degree = 3)
            : degree(degree), knots(knotsParam)
    {
        sort(begin(knots), end(knots));
        for(int i = 0; i < degree; ++i) {
            knots.insert(begin(knots), knots[0]);
            knots.push_back(knots.back());
        }
    }

    virtual ~BSpline() {
    }

    BSpline<T> copy() const override {
        return BSpline(getNonEndpointKnots(), degree);
    }

    vector<T> getNonEndpointKnots() const {
        vector<T> nonEndpointKnots;
        for(int i = degree; i < knots.size() - degree; ++i) {
            nonEndpointKnots.push_back(knots[i]);
        }
        return nonEndpointKnots;
    }

    vector<T> getCoordinates() const override {
        return getNonEndpointKnots();
    }

    void updateCoordinates(vector<T> const & direction, T stepLength) override {
        for(int i = degree; i < knots.size() - degree; ++i) {
            knots[i] -= stepLength * direction[i];
        }
    }

    BSpline<T> crossover(BSpline<T> const & g, float crossoverRate) const override {
        vector<T> knots1 = getNonEndpointKnots();
        vector<T> knots2 = g.getNonEndpointKnots();
        for(int i = 0; i < knots1.size(); ++i) {
            if(Util::rand::uniformReal() < crossoverRate) {
                knots1[i] = knots2[i];
            }
        }
        return BSpline(knots1, degree);
    }

    void mutate(float mutationRate, float relativeMutationSize) override {
        vector<T> nonEndpointKnots = getNonEndpointKnots();
        //Edge knots should stay fixed
        for(int i = 1; i < nonEndpointKnots.size()-1; ++i) {
            if(Util::rand::uniformReal() < mutationRate) {
                nonEndpointKnots[i] *= 1 + (2*Util::rand::uniformReal()-1) * relativeMutationSize;
            }
        }
        *this = BSpline(nonEndpointKnots, degree);
    }

    int getNumBases() const {
        return getNumBases(knots.size(), degree);
    }

    vector<T> getBasis(int i, vector<T> const & grid) const {
        auto basis = [&](T x){
            return BSpline::evaluateSpline(i, degree, knots, x);
        };
        vector<T> basisVector(grid.size());
        for(int j = 0; j < grid.size(); ++j) {
            basisVector[j] = grid[j] <= knots[i+degree+1] && grid[j] >= knots[i] ? basis(grid[j]) : 0;
        }
        return basisVector;
    }

    double evaluate(int i, T x) const {
        return evaluateSpline(i, degree, knots, x);
    }

    double mass(int basisIdx) const {
        int numPoints = 1000;
        T start = knots[basisIdx];
        T end = knots[basisIdx + degree + 1];
        T h = (end - start) / numPoints;
        T mass = 0;
        T x = start + h;
        for(int j = 0; j < numPoints; ++j) {
            mass += evaluate(basisIdx, x);
            x += h;
        }
        return mass * h;
    }

    vector<pair<T, T>> plotData(vector<T> const & grid, int basisIdx) const {
        vector<pair<T, T>> ret;
        ret.reserve(grid.size());
        for(T x : grid) {
            ret.emplace_back(x, evaluate(basisIdx, x));
        }
        return ret;
    }

    T max(int basisIdx) const {
        int numPoints = 1000;
        T start = knots[basisIdx];
        T end = knots[basisIdx + degree + 1];
        T h = (end - start) / numPoints;
        T maxValue = 0;
        T argMax = 0;
        T x = start + h;
        for(int j = 0; j < numPoints; ++j) {
            T y = evaluate(basisIdx, x);
            if (y > maxValue) {
                maxValue = y;
                argMax = x;
            }
            x += h;
        }
        return argMax;
    }

    void plotAllBases(string filenamePrefix, vector<T> const & grid) const {
        for(int basisIdx = 0; basisIdx < getNumBases(); ++basisIdx) {
            stringstream sstr;
            Util::writePlotData(plotData(grid, basisIdx), filenamePrefix, basisIdx);
        }
    }

    void updateInfo() const override {
        cout << "knots: ";
        for(T x : knots) cout << x << " ";
        cout << "\n";
        cout << "num bases: " << getNumBases() << "\n";
        for(int i = 0; i < getNumBases(); ++i) {
            cout << "max for basis " << i << ": " << max(i) << "\n";
        }
        vector<T> grid;
        vector<T> uniqueKnots;
        for(T x : knots) {

        }
        plotAllBases("spline", grid);
    }
};

#endif //CPP_BSPLINE_H
