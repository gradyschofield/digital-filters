//
// Created by Grady Schofield on 10/9/21.
//

#ifndef CPP_DISKRATIONAL_H
#define CPP_DISKRATIONAL_H

#include<vector>
#include<complex>
#include<cmath>

#include"Polynomial.h"
#include"Util.h"

using namespace std;

template<typename T>
class DiskRational {
    static const T constexpr maxNumeratorRadius = 1.2;
    static const T constexpr maxDenominatorRadius = 0.5;

    vector<T> numeratorCoef;
    vector<T> denominatorX;
    vector<T> denominatorTheta;

    vector<complex<T>> numeratorRoots;
    vector<complex<T>> denominatorRoots;
public:

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
        roots.reserve(denominatorTheta.size());
        for(int i = 0; i < denominatorTheta.size(); ++i) {
            roots.push_back(maxDenominatorRadius / (1 + exp(-denominatorX[i])) * exp(complex<T>(0, denominatorTheta[i])));
        }
        denominatorRoots = move(roots);
    }

    static DiskRational<T> create(int numNumeratorRoots, int numDenominatorRoots){
        Polynomial<T> numerator = Polynomial<T>::filterPolynomialFromRoots(numNumeratorRoots, maxNumeratorRadius);
        vector<complex<T>> denominatorRoots = Util::filterRoots(numDenominatorRoots, maxDenominatorRadius);
        return DiskRational(numerator.getCoefficients(), denominatorRoots);
    }

    DiskRational<T> copy() const {
        return DiskRational<T>(numeratorCoef, denominatorRoots);
    }

    void incorporateRoots(vector<complex<T>> numeratorRoots, vector<complex<T>> denominatorRoots) {
        Polynomial<T> p(numeratorCoef);
        p.incorporateRoots(numeratorRoots);
        numeratorCoef = p.getCoefficients();
        vector<T> newX;
        vector<T> newTheta;
        tie(newX, newTheta) = computeParamsFromRootsHelper(denominatorRoots, maxDenominatorRadius);
        for(T x : newX) denominatorX.push_back(x);
        for(T theta : newTheta) denominatorTheta.push_back(theta);
        computeRootsFromParams();
    }


    complex<T> evaluate(complex<T> x) const {
        complex<T> n = 0;
        complex<T> d = 1;
        complex<T> t = 1;
        for(T a : numeratorCoef) {
            n += a * t;
            t *= x;
        }
        for(complex<T> r : denominatorRoots) {
            d *= x - r;
            d *= x - conj(r);
        }
        return n / d;
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
        }
        computeRootsFromParams();
    }

    vector<T> getCoordinates() {
        vector<T> ret;
        std::copy(begin(numeratorCoef), end(numeratorCoef), back_inserter(ret));
        std::copy(begin(denominatorX), end(denominatorX), back_inserter(ret));
        std::copy(begin(denominatorTheta), end(denominatorTheta), back_inserter(ret));
        return ret;
    }

    tuple<vector<T>, T> derivative(vector<T> const & grid, vector<T> const & filter, int sampleRate) {
        vector<T> deriv(numeratorCoef.size() + denominatorX.size() + denominatorTheta.size());
        for(int i = 0; i < grid.size(); ++i) {
            T freq = grid[i];
            T filterValue = filter[i];
            complex<T> z = exp(complex<T>(0, -2 * M_PI * freq / sampleRate));
            complex<T> y = evaluate(z);
            complex<T> d = evaluateDenominator(z);
            T f = abs(y);
            T prefactor = 2 * (f - filterValue) / f;
            int k = 0;
            for(int i = 0; i < numeratorCoef.size(); ++i) {
                deriv[k] += prefactor * (conj(y) * Util::pow(z, i) / d).real();
                k += 1;
            }
            for(int i = 0; i < denominatorX.size(); ++i) {
                T x = denominatorX[i];
                T theta = denominatorTheta[i];
                complex<T> r = maxDenominatorRadius / (1 + exp(-x)) * exp(complex<T>(0, theta));
                complex<T> z2 = z - r;
                complex<T> z2rconj = z - conj(r);
                complex<T> dX = maxDenominatorRadius * exp(-x) / (T)pow(1 + exp(-x), 2) * exp(complex<T>(0, theta));
                complex<T> dTheta = maxDenominatorRadius * complex<T>(0, 1) * exp(complex<T>(0, theta)) / (1 + exp(-x));
                deriv[k] += prefactor * (conj(y) * y * (dX / z2 + conj(dX) / z2rconj)).real();
                deriv[k + 1] += prefactor * (conj(y) * y * (dTheta / z2 + conj(dTheta) / z2rconj)).real();
                k += 2;
            }
        }
        T numGridPointsInv = 1.0 / grid.size();
        T norm = 0;
        for(T & x : deriv) {
            x *= numGridPointsInv;
            norm += x*x;
        }
        return forward_as_tuple(move(deriv), sqrt(norm));
    }

    void mutate(float chanceOfMutation, float mutationSize) {
        for(T &x : numeratorCoef) {
            if(Util::rand::uniformReal() < chanceOfMutation) {
                x *= 1 + (Util::rand::uniformReal()-0.5)*mutationSize;
            }
        }
        for(T & x : denominatorX) {
            if(Util::rand::uniformReal() < chanceOfMutation) {
                x *= 1 + (Util::rand::uniformReal()-0.5)*mutationSize;
            }
        }
        for(T & theta : denominatorTheta) {
            if(Util::rand::uniformReal() < chanceOfMutation) {
                theta *= 1 + (Util::rand::uniformReal()-0.5)*mutationSize;
            }
        }
        computeRootsFromParams();
    }

    DiskRational crossover(DiskRational<T> const & other) const & {
        vector<T> newNumeratorCoef;
        newNumeratorCoef.reserve(numeratorCoef.size());
        for(int i = 0; i < numeratorCoef.size(); ++i) {
            newNumeratorCoef.push_back(Util::rand::uniformReal() < 0.5 ? numeratorCoef[i] : other.numeratorCoef[i]);
        }
        vector<complex<T>> newDenominatorRoots;
        newDenominatorRoots.reserve(denominatorRoots.size());
        for(int i = 0; i < denominatorRoots.size(); ++i) {
            newDenominatorRoots.push_back(Util::rand::uniformReal() < 0.5 ? denominatorRoots[i] : other.denominatorRoots[i]);
        }
        return DiskRational(newNumeratorCoef, newDenominatorRoots);
    }
};

#endif //CPP_DISKRATIONAL_H
