//
// Created by Grady Schofield on 10/17/21.
//

#include<Accelerate/Accelerate.h>

void floatBLASMultiply(int m, int n, int k, float alpha, float const * a, int lda, float const * b, int ldb, float beta, float * c, int ldc){
    cblas_sgemm(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans, CBLAS_TRANSPOSE::CblasNoTrans,
                m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

void doubleBLASMultiply(int m, int n, int k, double alpha, double const * a, int lda, double const * b, int ldb, double beta, double * c, int ldc){
    cblas_dgemm(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans, CBLAS_TRANSPOSE::CblasNoTrans,
                m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

