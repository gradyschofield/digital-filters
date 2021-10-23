//
// Created by Grady Schofield on 10/23/21.
//

#ifndef CPP_OPTIMIZABLE_H
#define CPP_OPTIMIZABLE_H

#include<vector>

using namespace std;

template<template<typename> typename R, typename T>
class Optimizable {
public:
    virtual R<T> copy() const = 0;
    virtual vector<T> getCoordinates() const;
    virtual void updateCoordinates(vector<T> const & direction, T stepLength) = 0;
    virtual ~Optimizable<R, T>(){
    }
};

#endif //CPP_OPTIMIZABLE_H
