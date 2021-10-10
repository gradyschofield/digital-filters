//
// Created by Grady Schofield on 10/9/21.
//

#ifndef CPP_POLYNOMIAL_H
#define CPP_POLYNOMIAL_H

#include<vector>
#include<complex>

#include<util.h>

using namespace std;

template<typename T>
class Polynomial {
    vector<T> coefficients;
public:
    Polynomial(vector<T> && coef)
        : coefficients(move(coef))
    {
    }

    /*
    Polynomial<T> mulMonomial(double a) {
        vector<T> newCoef = Util::shiftUp(coefficients);
        for(int i = 0; i < newCoef.size() - 1; ++i) {
            newCoef[i] += a * newCoef[i+1];
        }
        return Polynomial(move(newCoef));
    }

    static Polynomial<T> fromRoots(vector<T> const & roots) {
        Polynomial ret({1});
        for(T r : roots) {
            ret = ret.mulMonomial(-r)
        }
        return ret;
    }
     */

    Polynomial<T> incorporateRoots(vector<complex<T>> roots) {
        for(complex<T> r : roots) {
            vector<T> newCoef;
            newCoef.reserve(2 + coefficients.size());
            for (int i = 0; i < coefficients.size(); ++i) {
                newCoef[i + 2] = coefficients[i];
                newCoef[i + 1] -= 2 * r.real() * coefficients[i];
                newCoef[i] += norm(r) * coefficients[i];
            }
            coefficients = move(newCoef);
        }
    }

    vector<T> getCoefficients() const {
        return coefficients;
    }
};

#endif //CPP_POLYNOMIAL_H
