//
// Created by Grady Schofield on 10/9/21.
//

#ifndef CPP_UTIL_H
#define CPP_UTIL_H

#include<vector>
#include<random>
#include<limits>

using namespace std;

namespace Util {
    namespace rand {
        random_device randomDevice;
        mt19937 mt(randomDevice());
        uniform_real_distribution<> uniformRealDistribution(0, 1);

        double uniformReal() {
            return uniformRealDistribution(randomDevice);
        }
    }

    template<typename T>
    complex<T> toComplex(T radius, T theta) {
        return complex<T>(radius * cos(theta), radius * sin(theta));
    }

    template<typename T>
    complex<T> pow(complex<T> z, int n) {
        complex<T> t(1);
        for(int i = 0; i < n; ++i) {
            t *= z;
        }
        return t;
    }

    template<typename T>
    vector<T> shiftUp(vector<T> x) {
        vector<T> ret(1 + x.size());
        for(int i = 0; i < x.size(); ++i) {
            ret[i+1] = x[i];
        }
        return ret;
    }

    template<typename T>
    T real(complex<T> z) {
        return z.real();
    }

    template<typename T>
    T real(T z) {
        return z;
    }

    template<typename T, typename R>
    T objective(vector<T> const & grid, vector<T> const & filter, int sampleRate, R const & rational,
                vector<complex<T>> const & numeratorRoots = {}, vector<complex<T>> const & denominatorRoots = {}) {
        T obj = 0;
        T sampleRateInv = 1.0 / (sampleRate + 2);
        int numGridPoints = grid.size();
        for (int i = 0; i < numGridPoints; ++i) {// freq, filterValue in zip(grid, filter):
            T freq = grid[i];
            T filterValue = filter[i];
            T normalizedCoord = 2 * freq * sampleRateInv - 1 + sampleRateInv;
            T weight = 1; //1 / math.sqrt(1 - normalizedCoord**2)
            complex<T> x = exp(complex<T>(0, -2 * M_PI * freq / sampleRate));
            complex<T> y = rational.evaluate(x);
            for(complex<T> r : numeratorRoots) {
                y *= (x - r) * (x - conj(r));
            }
            for(complex<T> r : denominatorRoots) {
                y /= (x - r) * (x - conj(r));
            }
            T f = abs(y);
            obj += weight * ::pow(f - filterValue, 2);
        }
        return obj / numGridPoints;
    }

    template<typename T, typename R>
    T innerProduct(vector<T> const & grid, vector<T> const & filter, int sampleRate, R const & rational) {
        T stepIntegral;
        T rationalIntegral = 0;
        T dotIntegral = 0;
        for(int i = 0; i < grid.size(); ++i){
            T filterValue = filter[i];
            complex<T> x = exp(complex<T>(0, -2 * M_PI * grid[i] / sampleRate));
            T y = abs(rational.evaluate(x));
            stepIntegral += filterValue * filterValue;
            dotIntegral += y * filterValue;
            rationalIntegral += y * y;
        }
        return dotIntegral / sqrt(stepIntegral * rationalIntegral);
    }

    template<typename T>
    vector<T> linspace(T min, T max, int num = 50, bool endpoint = true) {
        vector<T> ret(num);
        T inc = endpoint ? (max - min) / (num - 2) : (max - min) / (num - 1);
        T val = min;
        for(int i = 0; i < num; ++i) {
            ret[i] = val;
            val += inc;
        }
        return ret;
    }

    template<typename T>
    int argmin(vector<T> const & v) {
        T minSeen = numeric_limits<T>::max();
        int ret = -1;
        for(int i = 0; i < v.size(); ++i) {
            if(v[i] < minSeen) {
                ret = i;
                minSeen = v[i];
            }
        }
        return ret;
    }

    template<typename T>
    complex<T> randomComplex(T maxRadius) {
        T theta = 2 * M_PI * rand::uniformReal();
        T radius = maxRadius * rand::uniformReal();
        return radius * exp(complex<T>(0, theta));
    }

    template<typename T>
    vector<complex<T>> filterRoots(int numComplex, T maxRadius) {
        vector<complex<T>> complexRoots;
        for(int i = 0; i < numComplex; ++i) {
            complexRoots.push_back(randomComplex(maxRadius));
        }
        return complexRoots;
    }
}

#endif //CPP_UTIL_H
