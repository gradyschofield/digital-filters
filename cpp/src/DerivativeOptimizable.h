//
// Created by Grady Schofield on 10/23/21.
//

#ifndef CPP_DERIVATIVEOPTIMIZABLE_H
#define CPP_DERIVATIVEOPTIMIZABLE_H

#include<vector>

using namespace std;

template<template<typename> typename Optimizable, typename FieldElement>
class DerivativeOptimizable {
public:
    virtual Optimizable<FieldElement> copy() const = 0;
    virtual vector<FieldElement> getCoordinates() const = 0;
    virtual void updateCoordinates(vector<FieldElement> const & direction, FieldElement stepLength) = 0;
    virtual ~DerivativeOptimizable<Optimizable, FieldElement>(){
    }
};

#endif //CPP_DERIVATIVEOPTIMIZABLE_H
