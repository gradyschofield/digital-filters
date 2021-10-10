//
// Created by Grady Schofield on 10/9/21.
//

#ifndef CPP_DISKRATIONAL_H
#define CPP_DISKRATIONAL_H

#include<vector>
#include<complex>
#include<cmath>

#include"Polynomial.h"

using namespace std;

template<typename T>
class DiskRational {
    static const T maxNumeratorRadius = 1.2;
    static const T maxDenominatorRadius = 0.5;

    vector<T> numeratorCoef;
    vector<T> denominatorX;
    vector<T> denominatorTheta;

    vector<complex<T>> numeratorRoots;
    vector<complex<T>> denominatorRoots;

    DiskRational(vector<T> const & numeratorCoef, vector<complex<T>> const & denominatorRoots)
        : numeratorCoef(numeratorCoef)
    {
        DiskRational::denominatorRoots = denominatorRoots;
        computeParamsFromRoots();
    }

    static tuple<vector<T>, vector<T>> computeParamsFromRootsHelper(vector<complex<T>> const & roots, T maxRadius) {
        vector<T> x;
        vector<T> theta;
        for(complex<T> r : roots) {
            x.push_back(-log(maxRadius / abs(r) - 1));
            theta.push_back(atan2(r.imag(), r.real()));
        }
        return forward_as_tuple(x, theta);
    }

    void computeParamsFromRoots() {
        tie(denominatorX, denominatorTheta) = computeParamsFromRootsHelper(denominatorRoots, maxDenominatorRadius);
    }

    /*
    def getRoots(self)

    :
    rc = list(self.numeratorCoef)
    rc.

    reverse()

    numeratorRoots = numpy.roots(rc)
    denominatorRoots =
    [
    complex(x,
    0) for
    x in
    self.realDenominatorRoots]
    denominatorRoots.
    extend(self
    .complexDenominatorRoots)
    return numeratorRoots,
    denominatorRoots

            def
    getAllRoots(self):
            rc = list(self.numeratorCoef)
    rc.

    reverse()

    numeratorRoots = numpy.roots(rc)
    denominatorRoots =
    [
    complex(x,
    0) for
    x in
    self.realDenominatorRoots]
    denominatorRoots.
    extend(self
    .complexDenominatorRoots)
    denominatorRoots.extend([z.

    conjugate()

    for
    z in
    self.complexDenominatorRoots])
    return numeratorRoots, denominatorRoots

     */

    void computeRootsFromParams() {
        vector<complex<T>> roots;
        roots.reserve(theta.size());
        for(int i = 0; i < theta.size(); ++i) {
            roots.push_back(maxDenominatorRadius / (1 + exp(-denominatorX[i])) * exp(complex(0, denominatorTheta[i]));
        }
        denominatorRoots = move(roots);
    }

    /*
    @
    staticmethod
            def
    create(numRealNumeratorRoots, numComplexNumeratorRoots, numRealDenominatorRoots, numComplexDenominatorRoots
    ):
    numerator = Util.filterPolynomialFromRoots(numRealNumeratorRoots, numComplexNumeratorRoots,
                                               DiskRational.maxNumeratorRadius)
    realDenominatorRoots,
    complexDenominatorRoots = Util.filterRoots(numRealDenominatorRoots, numComplexDenominatorRoots,
                                               DiskRational.maxDenominatorRadius)
    return

    DiskRational(numerator

    .coefficients, realDenominatorRoots, complexDenominatorRoots)

     */

    DiskRational copy() {
        return DiskRational(numeratorCoef, denominatorRoots);
    }

    void incorporateRoots(vector<complex<T>> numeratorRoots, vector<complex<T>> denominatorRoots) {
        Polynomial p(numeratorCoef);
        p.incorporateRoots(numeratorRoots);
        numeratorCoef = p.getCoefficients();
        vector<T> newX;
        vector<T> newTheta;
        tie(newX, newTheta) = computeParamsFromRootsHelper(denominatorRoots, maxDenominatorRadius);
        for(T x : newX) denominatorX.push_back(x);
        for(T theta : newTheta) denominatorTheta.push_back(theta);
        computeRootsFromParams()
    }


    complex<T> evaluate(complex<T> x) {
        complex<T> n = 0;
        complex<T> d = 1;
        T t = 1;
        for(T a : numeratorCoef) {
            n += a * t;
            t *= x;
        }
        for(complex<T> r : denominatorRoots) {
            d *= x - r;
            d *= x - conj(r);
        }
        return n / d
    }

    complex<T> evaluateDenominator(complex<T> x) {
        complex<T> d = 1;
        for(complex<T> r : denominatorRoots) {
            d *= x - r;
            d *= x - conj(r);
        }
        return d;
    }

    void updateCoef(vector<T> const & deriv, T rate) {
        int k = 0;
        for (T & c : numeratorCoef) {
            c -= rate * deriv[k];
            k += 1;
        }
        for (int i = 0; i < denominatorX.size(); ++i) {
            denominatorX[i] -= rate * deriv[k];
            denominatorTheta[i] -= rate * deriv[k + 1];
            T t = denominatorTheta[i] / (2 * M_PI);
            denominatorTheta[i] = (t - (int) (t)) * 2 * M_PI;
            k += 2;
            computeRootsFromParams();
        }
    }

    vector<T> getCoordinates() {
        vector<T> ret;
        copy(begin(numeratorCoef), end(numeratorCoef), back_inserter(ret));
        copy(begin(denominatorX), end(denominatorX), back_inserter(ret));
        copy(begin(denominatorTheta), end(denominatorTheta), back_inserter(ret));
        return ret;
    }

    vector<T> derivative(vector<T> const & grid, vector<T> const & filter, int sampleRate) {
        vector<T> deriv(numeratorCoef.size() + denominatorX.size() + denominatorTheta.size());
        for(int i = 0; i < grid.size(); ++i) {
            T freq = grid[i];
            T filterValue = filter[i];
            complex<T> z = exp(complex(0, -2 * math.pi * freq / sampleRate));
            complex<T> y = evaluate(z);
            complex<T> d = evaluateDenominator(z);
            T f = abs(y);
            T prefactor = 2 * (f - filterValue) / f;
            k = 0
            for(int i = 0; i < numeratorCoef.size(); ++i) {
                deriv[k] += prefactor * (conj(y) * fixedPart * pow(z, i) / d).real()
                k += 1
            }
            for i in range(len(self .realDenominatorX)):
                x = self.realDenominatorX[i]
                r = DiskRational.maxDenominatorRadius * (1 - math.exp(-x) / (1 + math.exp(-x)))
                z2 = z - r
                deriv[k] += prefactor * DiskRational . maxDenominatorRadius * (y * y / z2 * fixedPart * (-2 * math.exp(-2 * x) / (1 + math.exp(-x)) * *2)). real
                k += 1
            for x, theta in zip(self .complexDenominatorX, self.complexDenominatorTheta):
            r = DiskRational.maxDenominatorRadius / (1 + math.exp(-x)) * cmath.exp(complex(0, theta))
            z2 = z - r
            z2rconj = z - r.conjugate()
            dX = DiskRational.maxDenominatorRadius * math.exp(-x) / (1 + math.exp(-x)) * *2 *
                 cmath.exp(complex(0, theta))
            dTheta = DiskRational.maxDenominatorRadius * complex(0, 1) * cmath.exp(complex(0, theta)) /
                     (1 + math.exp(-x))
            deriv[k] +=
                    prefactor * (y
                                         .

                                                 conjugate()

                                 *
                                 fixedPart * y
                                 * (dX / z2 + dX.

                            conjugate()

                                              / z2rconj)).
                            real
            deriv[k + 1]
                    +=
                    prefactor * (y
                                         .

                                                 conjugate()

                                 *
                                 fixedPart * y
                                 * (dTheta / z2 + dTheta.

                            conjugate()

                                                  / z2rconj)).
                            real
            k
                    += 2
        }
        numGridPoints = len(grid)
        deriv =
        [t / numGridPoints
        for
        t in
        deriv]
        norm = functools.reduce(lambda
        n,
                t: n + t * t, deriv,
                0)
        return deriv, math.
                sqrt(norm)
    }

    def mutate(self, chanceOfMutation, mutationSize)

    :
    newCoef = []
    for
    x in
    self.
    numeratorCoef:
    if random.uniform(0, 1) <
    chanceOfMutation:
            x
    *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
    newCoef.
    append(x)
    self.
    numeratorCoef = newCoef
    newCoef = []
    for
    x in
    self.
    realDenominatorX:
    if random.uniform(0, 1) <
    chanceOfMutation:
            x
    *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
    newCoef.
    append(x)
    self.
    realDenominatorX = newCoef
    newCoef = []
    for
    x in
    self.
    complexDenominatorX:
    if random.uniform(0, 1) <
    chanceOfMutation:
            x
    *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
    newCoef.
    append(x)
    self.
    complexDenominatorX = newCoef
    newCoef = []
    for
    x in
    self.
    complexDenominatorTheta:
    if random.uniform(0, 1) <
    chanceOfMutation:
            x
    *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
    newCoef.
    append(x)
    self.
    complexDenominatorTheta = newCoef
    self.

    computeRootsFromParams()


    def crossover(self, other)

    :
    newNumeratorCoef = []
    for x,
    y in
    zip(self
    .numeratorCoef, other.numeratorCoef):
    newNumeratorCoef.
    append(x
    if random.uniform(0, 1) < 0.5 else y)
    newRealDenominatorRoots = []
    for x,
    y in
    zip(self
    .realDenominatorRoots, other.realDenominatorRoots):
    newRealDenominatorRoots.
    append(x
    if random.uniform(0, 1) < 0.5 else y)
    newComplexDenominatorRoots = []
    for x,
    y in
    zip(self
    .complexDenominatorRoots, other.complexDenominatorRoots):
    newComplexDenominatorRoots.
    append(x
    if random.uniform(0, 1) < 0.5 else y)
    return

    DiskRational(newNumeratorCoef, newRealDenominatorRoots, newComplexDenominatorRoots)
};

#endif //CPP_DISKRATIONAL_H
