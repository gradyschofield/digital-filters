//
// Created by Grady Schofield on 10/17/21.
//

#ifndef CPP_BLAS_H
#define CPP_BLAS_H

template<typename T>
void blasMultiply(int m, int n, int k, T alpha, T const * a, int lda, T const * b, int ldb, T beta, T * c, int ldc);

template<typename T>
vector<complex<T>> blasEigenvalues(int n, T const * a, int lda);

void floatBLASMultiply(int m, int n, int k, float alpha, float const * a, int lda, float const * b, int ldb, float beta, float * c, int ldc);
void doubleBLASMultiply(int m, int n, int k, double alpha, double const * a, int lda, double const * b, int ldb, double beta, double * c, int ldc);
tuple<vector<float>, vector<float>> floatEigenvalues(int n, float * a, int lda);
tuple<vector<double>, vector<double>> doubleEigenvalues(int n, double * a, int lda);

template<>
void blasMultiply(int m, int n, int k, float alpha, float const * a, int lda, float const * b, int ldb, float beta, float * c, int ldc) {
    floatBLASMultiply(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

template<>
void blasMultiply(int m, int n, int k, double alpha, double const * a, int lda, double const * b, int ldb, double beta, double * c, int ldc){
    doubleBLASMultiply(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

#include<Complex.h>

template<>
vector<complex<float>> blasEigenvalues(int n, float const * a, int lda) {
    vector<float> re;
    vector<float> im;
    vector<float> aCopy(n*lda);
    memcpy(aCopy.data(), a, n*lda*sizeof(float));
    tie(re, im) = floatEigenvalues(n, aCopy.data(), lda);
    vector<complex<float>> ret;
    for(int i = 0; i < n; ++i) {
        ret.emplace_back(re[i], im[i]);
    }
    return ret;
}

template<>
vector<complex<double>> blasEigenvalues(int n, double const * a, int lda) {
    vector<double> re;
    vector<double> im;
    vector<double> aCopy(n*lda);
    memcpy(aCopy.data(), a, n*lda*sizeof(double));
    tie(re, im) = doubleEigenvalues(n, aCopy.data(), lda);
    vector<complex<double>> ret;
    for(int i = 0; i < n; ++i) {
        ret.emplace_back(re[i], im[i]);
    }
    return ret;
}

#endif //CPP_BLAS_H
