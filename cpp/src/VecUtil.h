//
// Created by Grady Schofield on 10/10/21.
//

#ifndef CPP_VECUTIL_H
#define CPP_VECUTIL_H


#include<stdexcept>
#include<typeinfo>

#include<BLAS.h>

using namespace std;

namespace VecUtil {
    template<typename T>
    T norm(vector<T> const & v) {
        T n = 0;
        for(T x : v) {
            n += x*x;
        }
        return sqrt(n);
    }

    /*
    template<typename T>
    T dot(vector<complex<T>> const & v1, vector<complex<T>> const &v2) {
        T d = 0;
        for(int i = 0; i < v1.size(); ++i) {
            d += conj(v1[i]) * v2[i];
        }
        return d;
    }
     */

    template<typename T>
    T dot(vector<T> const & v1, vector<T> const &v2) {
        T d = 0;
        for(int i = 0; i < v1.size(); ++i) {
            d += v1[i] * v2[i];
        }
        return d;
    }

    template<typename T>
    vector<T> subtract(vector<T> const & v1, vector<T> const &v2) {
        vector<T> ret(v1);
        for(int i = 0; i < v1.size(); ++i) {
            ret[i] -= v2[i];
        }
        return ret;
    }

    template<typename T>
    class Vector {
        Vector<T> operator-(Vector<T> const &v) {
        }

        T dot(Vector<T> const &v) {

        }
    };

    template<typename T>
    class Matrix {
        int numRows, numCols, leadingDimension;
        vector<T> matrix;

        int index(int i, int j) const {
            return j + i*numCols;
        }

    public:
        Matrix(int numRows, int numCols)
                : numRows(numRows), numCols(numCols), leadingDimension(numRows), matrix(numRows * numCols)
        {
        }

        Matrix(int numRows)
                : numRows(numRows), numCols(numRows), leadingDimension(numRows), matrix(numRows * numCols)
        {
        }

        int getNumRows() const {
            return numRows;
        }

        int getNumCols() const {
            return numCols;
        }

        template<typename Func>
        void fill(Func && fillFunc) {
            for(int i = 0; i < numRows; ++i) {
                for(int j = 0; j < numCols; ++j) {
                    matrix[index(i, j)] = fillFunc(i, j);
                }
            }
        }

        Matrix<T> transpose() const {
            Matrix<T> ret(numRows, numCols);
            for(int i = 0; i < numRows; ++i) {
                for (int j = 0; j < numCols; ++j) {
                    ret.matrix[index(j,i)] = matrix[index(i, j)];
                }
            }
            return ret;
        }

        vector<complex<T>> eigenvalues() const {
            return blasEigenvalues(numRows, matrix.data(), leadingDimension);
        }

        Matrix<T> multiply(Matrix<T> const & m) const {
            Matrix<T> ret(numRows, numCols);
            if(true) {
                blasMultiply(numRows, m.numCols, numCols, (T)1, matrix.data(), numRows, m.matrix.data(), m.numRows,
                             (T)0, ret.matrix.data(), numRows);
            } else {
                for (int i = 0; i < numRows; ++i) {
                    for (int j = 0; j < numCols; ++j) {
                        for (int k = 0; k < numCols; ++k) {
                            ret.matrix[index(i, j)] += matrix[index(i, k)] * m.matrix[index(k, j)];
                        }
                    }
                }
            }
            return ret;
        }

        vector<T> multiply(vector<T> const & v) const {
            vector<T> ret(numRows);
            for (int i = 0; i < numRows; ++i) {
                T d = 0;
                for (int j = 0; j < numCols; ++j) {
                    d += matrix[index(i, j)] * v[j];
                }
                ret[i] = d;
            }
            return ret;
        }

        Matrix<T> subtract(Matrix<T> const & m) const {
            Matrix<T> ret(numRows, numCols);
            ret.fill([this, &m](int i, int j) {
                return matrix[index(i,j)] - m.matrix[index(i,j)];
            });
            return ret;
        }

        Matrix<T> add(Matrix<T> const & m) const {
            Matrix<T> ret(numRows, numCols);
            ret.fill([this, &m](int i, int j) {
                return matrix[index(i,j)] - m.matrix[index(i,j)];
            });
            return ret;
        }

        Matrix<T> operator*(Matrix<T> const & m) {
            return multiply(m);
        }

        vector<T> operator*(vector<T> const & v) {
            return multiply(v);
        }

        Matrix<T> operator-(Matrix<T> const & m) {
            return subtract(m);
        }

        Matrix<T> operator+(Matrix<T> const & m) {
            return add(m);
        }

        static Matrix<T> identity(int numRows) {
            Matrix ret(numRows, numRows);
            ret.fill([](int i, int j) {
                return i == j ? T(1) : T(0);
            });
            return ret;
        }
    };
}

#endif //CPP_VECUTIL_H
