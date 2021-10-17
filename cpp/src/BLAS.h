//
// Created by Grady Schofield on 10/17/21.
//

#ifndef CPP_BLAS_H
#define CPP_BLAS_H

template<typename T>
void blasMultiply(int m, int n, int k, T alpha, T const * a, int lda, T const * b, int ldb, T beta, T * c, int ldc);

void floatBLASMultiply(int m, int n, int k, float alpha, float const * a, int lda, float const * b, int ldb, float beta, float * c, int ldc);
void doubleBLASMultiply(int m, int n, int k, double alpha, double const * a, int lda, double const * b, int ldb, double beta, double * c, int ldc);

template<>
void blasMultiply(int m, int n, int k, float alpha, float const * a, int lda, float const * b, int ldb, float beta, float * c, int ldc) {
    floatBLASMultiply(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

template<>
void blasMultiply(int m, int n, int k, double alpha, double const * a, int lda, double const * b, int ldb, double beta, double * c, int ldc){
    doubleBLASMultiply(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

#endif //CPP_BLAS_H
