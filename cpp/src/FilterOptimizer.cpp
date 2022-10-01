/*
 *
 * This program finds an infinite impulse response filter whose magnitude response, |H(z)|, is optimized to match
 * a B-spline basis function.  The point of using B-splines is that they form a partition of unity.  The phase
 * response is ignored.
 *
 * ---------------------------------------------------------------------------------------------------------------------
 *
 * The Setup
 * =========
 * Splines are not used to represent the filter during the optimization.
 * Splines were chosen merely for their partition of unity property.  The filter function is represented
 * by its transfer function's roots in the complex plane.  These roots are the degrees of freedom (DOF) in the
 * optimization.  The objective function is evaluated by comparing the magnitude response to the ideal filter.
 * This means we are taking the complex roots of a rational polynomial as our DOF, but evaluating
 * the objective function on a grid in (real) frequency space.  When we refer to the "numerator" or "denominator",
 * we are talking about the transfer function H(z), specifically for a filter with outputs y_i and inputs x_i,
 *
 *      y_n = b_0 * x_n + b_1 * x_(n-1) + ... b_k * x_(n-k)
 *              - a_1 * y_(n-1) - a_2 * y_(n-2) - ... - a_j * y_(n-j)
 *
 * we have the transfer function,
 *
 *      H(z) = (b_0 * z + b_1 * z^-1 + ... + b_k * z^-k)
 *              / (1 + a_1 * y^-1 + ... + a_j * z^-j),
 *
 * where, for a given frequency f and sample rate (frequency) Fs,
 *
 *      z = exp(2i * pi * f / Fs).
 *
 * Note that the coefficients a_i and b_i are always real, so when we insert roots in the transfer function, they
 * come in conjugate pairs.
 *
 * If you want more details see Chapter 1 of Rusty Allred's book, Digital Filters for Everyone.
 *
 * ---------------------------------------------------------------------------------------------------------------------
 *
 * Algorithm
 * =========
 *
 * 1. Repeatedly reduce the uniform grid's spacing until it represents the spline basis function reasonably well.
 *    These grid points span the frequency range from 0 to half the sample rate Fs.
 * 2. Do a few steps of a genetic optimization routine to get good starting DOF.
 * 3. Use a polar grid in the complex plane to roughly find the best place to insert a root.  Separately try a numerator
 *    root and a denominator root.  Insert either the numerator or denominator but not both.
 * 4. Do a line search with the new DOF to refine the position of the roots
 * 5. Go back to step 3
 *
 *
 */
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
        vector<T> interpolatedFunc = Util::interpolate(grid, filterFunc, fineGrid);
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
    rational = lineSearch(objective, derivative, rational, 1E-5, 10, 10, 1E-7, true);
    while(true) {
        complex<T> minNumeratorRoot, minDenominatorRoot;
        T minNumeratorObj, minDenominatorObj;
        tie(minNumeratorRoot, minNumeratorObj) = Util::minimizeOnGrid(grid, filterFunc, sampleRate, rational, true);
        tie(minDenominatorRoot, minDenominatorObj) = Util::minimizeOnGrid(grid, filterFunc, sampleRate, rational, false);
        if (minNumeratorObj < minDenominatorObj) {
            rational.incorporateRoots({minNumeratorRoot}, {});
        } else {
            rational.incorporateRoots({}, {minDenominatorRoot});
        }
        rational = lineSearch(objective, derivative, rational, 1E-5, 10, 10, 1E-7, true);
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
