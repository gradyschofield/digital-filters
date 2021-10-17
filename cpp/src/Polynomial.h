//
// Created by Grady Schofield on 10/9/21.
//

#ifndef CPP_POLYNOMIAL_H
#define CPP_POLYNOMIAL_H

#include<vector>
#include<complex>

#include<Util.h>

using namespace std;

template<typename T>
class Polynomial {
    vector<T> coefficients;
public:
    Polynomial(vector<T> && coef)
        : coefficients(move(coef))
    {
    }

    Polynomial(vector<T> const & coef)
        : coefficients(coef)
    {
    }

    Polynomial<T> mulMonomial(T a) {
        vector<T> newCoef = Util::shiftUp(coefficients);
        for(int i = 0; i < newCoef.size() - 1; ++i) {
            newCoef[i] += a * newCoef[i+1];
        }
        return Polynomial<T>(move(newCoef));
    }

    static Polynomial<T> fromRoots(vector<T> const & roots) {
        Polynomial<T> ret({1});
        for(T r : roots) {
            ret = ret.mulMonomial(-r);
        }
        return ret;
    }

    void incorporateRoots(vector<complex<T>> roots) {
        for(complex<T> r : roots) {
            vector<T> newCoef(2 + coefficients.size());
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

    static Polynomial<T> filterPolynomialFromRoots(int numComplex, T radius) {
        vector<complex<T>> roots;
        for(int i = 0; i < numComplex; ++i) {
            complex<T> r = Util::randomComplex(radius);
            roots.push_back(r);
            roots.push_back(conj(r));
        }
        return Polynomial<complex<T>>::fromRoots(roots).template realifyCoef<T>();
    }

    template<typename R>
    Polynomial<R> realifyCoef() const {
        vector<R> newCoef;
        for(T z : coefficients) {
            newCoef.push_back(Util::real(z));
        }
        return Polynomial<R>(move(newCoef));
    }

};

#endif //CPP_POLYNOMIAL_H
