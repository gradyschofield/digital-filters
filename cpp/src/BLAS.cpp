//
// Created by Grady Schofield on 10/17/21.
//

#include<vector>
#include<Accelerate/Accelerate.h>

using namespace std;

void floatBLASMultiply(int m, int n, int k, float alpha, float const * a, int lda, float const * b, int ldb, float beta, float * c, int ldc){
    cblas_sgemm(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans, CBLAS_TRANSPOSE::CblasNoTrans,
                m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

void doubleBLASMultiply(int m, int n, int k, double alpha, double const * a, int lda, double const * b, int ldb, double beta, double * c, int ldc){
    cblas_dgemm(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans, CBLAS_TRANSPOSE::CblasNoTrans,
                m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

tuple<vector<float>, vector<float>> floatEigenvalues(int n, float * a, int lda){
    char job = 'N';
    vector<float> wr(n), wi(n);
    int ldvl = lda;
    vector<float> work(n*n);
    int workSize = work.size();
    int info;
    sgeev_(&job, &job, &n, a, &lda, wr.data(), wi.data(), nullptr, &ldvl, nullptr, &ldvl, work.data(), &workSize, &info);
    if(info != 0) {
        throw runtime_error("sgeev_ failed");
    }
    return forward_as_tuple(move(wr), move(wi));
}

tuple<vector<double>, vector<double>> doubleEigenvalues(int n, double * a, int lda){
    char job = 'N';
    vector<double> wr(n), wi(n);
    int ldvl = lda;
    vector<double> work(n*n);
    int workSize = work.size();
    int info;
    dgeev_(&job, &job, &n, a, &lda, wr.data(), wi.data(), nullptr, &ldvl, nullptr, &ldvl, work.data(), &workSize, &info);
    if(info != 0) {
        throw runtime_error("dgeev_ failed");
    }
    return forward_as_tuple(move(wr), move(wi));
}

