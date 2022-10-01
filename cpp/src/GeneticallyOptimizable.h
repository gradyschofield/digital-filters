//
// Created by Grady Schofield on 10/24/21.
//

#ifndef CPP_GENETICALLYOPTIMIZABLE_H
#define CPP_GENETICALLYOPTIMIZABLE_H

template<template<typename> typename Optimizable, typename FieldElement>
class GeneticallyOptimizable {
public:
    virtual Optimizable<FieldElement> crossover(Optimizable<FieldElement> const & g, float crossoverRate) const = 0;
    virtual void mutate(float mutationRate, float relativeMutationSize) = 0;
    virtual void updateInfo() const {
    }
    ~GeneticallyOptimizable(){
    }
};

#endif //CPP_GENETICALLYOPTIMIZABLE_H
