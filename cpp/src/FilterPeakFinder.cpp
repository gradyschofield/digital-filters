//
// Created by Grady Schofield on 10/23/21.
//

#include<iostream>

#include<BSpline.h>

int main(int argc, char ** argv) {

    typedef double T;

    int sampleRate = 96000;
    int numGridPoints = 10000;
    int basisIdx = 7;
    int degree = 3;

    T minFrequency = 0;
    T maxFrequency = sampleRate / 2;
    //vector<T> basisMaxima{39, 78, 156, 312, 625, 1250, 2500, 5000, 10000, 20000};
    vector<T> basisMaxima{0, 12.5, 25, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 25600*2, 25600*4,25600*8 };
    vector<T> maximumCeneteredBasisKnots{minFrequency};
    for(int i = 0; i < degree; ++i) {
        //maximumCeneteredBasisKnots.push_back(minFrequency + (i+1) * 0.01);
        //maximumCeneteredBasisKnots.push_back(minFrequency);
    }
    copy(begin(basisMaxima), end(basisMaxima), back_inserter(maximumCeneteredBasisKnots));
    for(int i = 0; i < degree; ++i) {
        T h = (maxFrequency - basisMaxima.back()) * 0.01;
        //maximumCeneteredBasisKnots.push_back(maxFrequency - (i+1) * 0.01);
        //maximumCeneteredBasisKnots.push_back(maxFrequency);
    }

    vector<T> grid = Util::linspace<T>(0, sampleRate / 2, numGridPoints);
    //BSpline<T> splineBasis(grid, BSpline<T>::getDefaultKnots(), degree);
    BSpline<T> splineBasis(grid, maximumCeneteredBasisKnots, degree);
    cout << "default knots: ";
    for(T x : BSpline<T>::getDefaultKnots()) cout << x << " ";
    cout << "\n";
    cout << "default knots: ";
    for(T x : maximumCeneteredBasisKnots) cout << x << " ";
    cout << "\n";
    for(int basisIdx = 0; basisIdx < splineBasis.getNumBases(); ++basisIdx) {
        stringstream sstr;
        Util::writePlotData(splineBasis.plotData(grid, basisIdx), "spline", basisIdx);
        cout << splineBasis.max(basisIdx) << "\n";
    }

    return 0;
}