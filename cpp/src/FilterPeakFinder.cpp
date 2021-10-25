//
// Created by Grady Schofield on 10/23/21.
//

#include<iostream>

#include<BSpline.h>
#include<GeneticAlgorithm.h>

int main(int argc, char ** argv) {

    typedef double T;

    int sampleRate = 96000;
    int numGridPoints = 10000;
    int degree = 3;

    vector<T> basisMaximum{5, 50, 75, 150, 300, 600, 900, 1250, 1800, 2500, 3500, 5000, 7500, 10000, 20000, 30000};
    //vector<T> knots{0, 5, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 40000, 48000};
    vector<T> knots{0, 5, 50, 75, 150, 300, 600, 900, 1250, 1800, 2500, 3500, 5000, 7500, 10000, 20000, 25000, 30000, 35000, 48000};

    auto create = [&knots, degree]() {
        BSpline<T> ret(knots, degree);
        ret.mutate(0.5, 0.01);
        return ret;
    };

    auto objective = [&basisMaximum](BSpline<T> const & spline) {
        int numBases = spline.getNumBases();
        T obj = 0;
        for(int i = 1; i < min(numBases, (int)basisMaximum.size()); ++i) {
            obj += fabs(spline.max(i) - basisMaximum[i-1])/basisMaximum[i-1];
        }
        return obj/min(numBases, (int)basisMaximum.size());
    };

    vector<T> grid = Util::linspace<T>(0, sampleRate / 2, numGridPoints);
    BSpline<T> splineBasis(knots, degree);
    splineBasis.printInfo();
    splineBasis.plotAllBases("spline", grid);
    BSpline<T> optimizedSpline = geneticOptimizer<BSpline, T>(objective, create, 10000);
    optimizedSpline.plotAllBases("spline", grid);
    optimizedSpline.printInfo();
    return 0;
}