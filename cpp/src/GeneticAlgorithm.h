//
// Created by Grady Schofield on 10/23/21.
//

#ifndef CPP_GENETICALGORITHM_H
#define CPP_GENETICALGORITHM_H

#include<vector>
#include<iostream>
#include<concepts>
#include<random>

#include<Util.h>
#include<GeneticallyOptimizable.h>

using namespace std;

template<template<typename> typename R, typename T, typename Objective, typename Create>
requires(is_floating_point_v<T> &&
         is_invocable_r_v<T, Objective, R<T> const &> &&
         is_invocable_r_v<R<T>, Create> &&
         is_base_of<GeneticallyOptimizable<R, T>, R<T>>::value)
R<T> geneticOptimizer(Objective && objective, Create && create,
                      int maxGenerations = 5, int populationSize = 2000, int cullSize = 200,
                      float crossoverRate = 0.5, float mutationRate = 0.5, float relativeMutationSize = 0.01) {
    vector<pair<R<T>, T>> population;
    for(int i = 0; i < populationSize - population.size(); ++i) {
        R<T> r = invoke(forward<Create>(create));
        T obj = invoke(forward<Objective>(objective), forward<R<T> const &>(r));
        population.emplace_back(move(r), obj);
    }

    sort(begin(population), end(population), [](auto & x, auto & y){
        return x.second < y.second;
    });
    population.erase(begin(population) + cullSize, end(population));
    cout << "best obj " <<  population[0].second << "\n";
    uniform_int_distribution<> uniformIntDistribution(0, cullSize-1);
    for(int j = 0; j < maxGenerations; ++j) {
        for(int i = cullSize; i < populationSize; ++i) {
            R<T> const & r1 = population[uniformIntDistribution(Util::rand::randomDevice)].first;
            R<T> const & r2 = population[uniformIntDistribution(Util::rand::randomDevice)].first;
            R<T> r = r1.crossover(r2, crossoverRate);
            r.mutate(mutationRate, relativeMutationSize);
            T obj = invoke(forward<Objective>(objective), forward<R<T> const &>(r));
            population.emplace_back(move(r), obj);
        }
        sort(begin(population), end(population), [](auto & x, auto & y){
            return x.second < y.second;
        });
        population.erase(begin(population) + cullSize, end(population));
        cout << "best obj " << population[0].second << "\n";
        if(j % 50 == 0) {
            population[0].first.updateInfo();
        }
    }
    T obj = invoke(forward<Objective>(objective), forward<R<T> const &>(population[0].first));
    cout << "obj: " << obj << "\n";
    return move(population[0].first);
}

#endif //CPP_GENETICALGORITHM_H
