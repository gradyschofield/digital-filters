//
// Created by Grady Schofield on 10/23/21.
//

#ifndef CPP_DERIVATIVEOPTIMIZABLE_H
#define CPP_DERIVATIVEOPTIMIZABLE_H

#include<vector>

using namespace std;

template<template<typename> typename R, typename T>
class DerivativeOptimizable {
public:
    virtual R<T> copy() const = 0;
    virtual vector<T> getCoordinates() const = 0;
    virtual void updateCoordinates(vector<T> const & direction, T stepLength) = 0;
    virtual ~DerivativeOptimizable<R, T>(){
    }
};

#endif //CPP_DERIVATIVEOPTIMIZABLE_H
