#include <iostream>
#include <vector>
#include <chrono>

#include<DiskRational.h>
#include<Util.h>
#include<VecUtil.h>
#include<BSpline.h>
#include<LineSearch.h>
#include<GeneticAlgorithm.h>

using namespace std;

template<typename T>
tuple<complex<T>, T> minimizeOnGrid(vector<T> const & grid, vector<T> const & filterFunc,
                                    int sampleRate, DiskRational<T> const & rational, bool insertNumerator) {
    T minObj = numeric_limits<T>::max();
    complex<T> minRoot;
    for (T theta: Util::linspace<T>(0, M_PI, 20)) {
        for (T radius: Util::linspace<T>(0.01, insertNumerator ? 1.1 : 0.5, 20)) {
            complex<T> root = Util::toComplex(radius, theta);
            T obj;
            if (insertNumerator) {
                obj = Util::objective(grid, filterFunc, sampleRate, rational, {root});
            } else {
                obj = Util::objective(grid, filterFunc, sampleRate, rational, {}, {root});
            }
            if (obj < minObj) {
                minObj = obj;
                minRoot = root;
            }
        }
    }
    return forward_as_tuple(minRoot, minObj);
}

template<typename T>
vector<T> interpolate(vector<T> const & coarseGrid, vector<T> const & coarseFunc, vector<T> const & fineGrid) {
    vector<T> fineFunc(fineGrid.size());
    int coarseGridIdx = 0;
    for(int i = 0; i < fineGrid.size(); ++i) {
        T x = fineGrid[i];
        while(x > coarseGrid[coarseGridIdx+1]) {
            ++coarseGridIdx;
        }
        T alpha = (x - coarseGrid[coarseGridIdx]) / (coarseGrid[coarseGridIdx+1] - coarseGrid[coarseGridIdx]);
        fineFunc[i] = (1-alpha) * coarseFunc[coarseGridIdx] + alpha * coarseFunc[coarseGridIdx+1];
    }
    return fineFunc;
}

int main() {

    typedef double T;

    int sampleRate = 96000;
    int numGridPoints = 50;
    int numNumeratorRoots = 20;
    int numDenominatorRoots = 20;
    int basisIdx = 8;

    bool foundGrid = false;
    for(; numGridPoints < 10000; numGridPoints += 50) {
        vector<T> grid = Util::linspace<T>(0, sampleRate/2, numGridPoints);
        BSpline<T> splineBasis(BSpline<T>::getDefaultKnots());
        vector<T> filterFunc = splineBasis.getBasis(basisIdx, grid);
        vector<T> fineGrid = Util::linspace<T>(0, sampleRate/2, numGridPoints + 50);
        BSpline<T> fineSplineBasis(BSpline<T>::getDefaultKnots());
        vector<T> fineFilterFunc = fineSplineBasis.getBasis(basisIdx, fineGrid);
        vector<T> interpolatedFunc = interpolate(grid, filterFunc, fineGrid);
        T basisMass = splineBasis.mass(basisIdx);
        T n = 0;
        T approximateMass = 0;
        T fineGridSpacing = fineGrid[1] - fineGrid[0];
        for(int j = 0; j < fineGrid.size(); ++j) {
            n += pow(splineBasis.evaluate(basisIdx, fineGrid[j]) - interpolatedFunc[j], 2);
            approximateMass += interpolatedFunc[j];
        }
        approximateMass *= fineGridSpacing;
        n = sqrt(n / fineGrid.size());
        cout << "fine grid " << numGridPoints << " norm: " << n <<
             " basis mass: " << basisMass << " approx mass: " << approximateMass << "\n";
        if(n < 5E-4 && fabs(basisMass-approximateMass)/basisMass < 0.02) {
            foundGrid = true;
            break;
        }
    }
    if(!foundGrid) {
        cout << "Couldn't find a reasonable grid spacing that produced a reasonably converged approximation to the filter function.\n";
        exit(1);
    }
    vector<T> grid = Util::linspace<T>(0, sampleRate/2, numGridPoints);
    BSpline<T> splineBasis(BSpline<T>::getDefaultKnots());
    vector<T> filterFunc = splineBasis.getBasis(basisIdx, grid);

    auto objective = [&grid, &filterFunc, sampleRate](DiskRational<T> const & r) {
        return Util::objective(grid, filterFunc, sampleRate, r);
    };

    auto derivative = [&grid, &filterFunc, sampleRate](DiskRational<T> const & r) {
        return r.derivative(grid, filterFunc, sampleRate);
    };

    auto create = [numNumeratorRoots, numDenominatorRoots]() {
        return DiskRational<T>::create(numNumeratorRoots, numDenominatorRoots);
    };

    DiskRational<T> rational = geneticOptimizer<DiskRational, T>(objective, create);
    ofstream rootTalk("rootTalk");
    auto startTime = chrono::high_resolution_clock::now();
    int iterations = 0;
    while(true) {
        rational = lineSearch(objective, derivative, rational, 1E-5, 10, 10, 1E-7, true);
        complex<T> minNumeratorRoot, minDenominatorRoot;
        T minNumeratorObj, minDenominatorObj;
        tie(minNumeratorRoot, minNumeratorObj) = minimizeOnGrid(grid, filterFunc, sampleRate, rational, true);
        tie(minDenominatorRoot, minDenominatorObj) = minimizeOnGrid(grid, filterFunc, sampleRate, rational, false);
        if (minNumeratorObj < minDenominatorObj) {
            rational.incorporateRoots({minNumeratorRoot}, {});
        } else {
            rational.incorporateRoots({}, {minDenominatorRoot});
        }
        ++iterations;
        Util::writePlotData(rational.plotData(grid, sampleRate), "approx", iterations);
        Util::writePlotData(rational.plotResidualData(grid, filterFunc, sampleRate), "residual", iterations);
        Util::writePlotData(rational.getNumeratorRoots(), "numeratorRoots", iterations);
        Util::writePlotData(rational.getDenominatorRoots(), "denominatorRoots", iterations);
        /*
        auto endTime = chrono::high_resolution_clock::now();
        double iterationsPerSecond = numDerivativesCalculated / (double)chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count();
        rootTalk << min(minNumeratorObj, minDenominatorObj) << " iterations per second " << 1E9*iterationsPerSecond <<
            " inserting in " << (minNumeratorObj < minDenominatorObj ? "numerator" : "denominator") << endl;
            */
    }

    return 0;
}
