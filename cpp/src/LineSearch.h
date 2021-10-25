//
// Created by Grady Schofield on 10/23/21.
//

#ifndef CPP_LINESEARCH_H
#define CPP_LINESEARCH_H

#include<concepts>
#include<vector>
#include<iostream>

#include<VecUtil.h>

using namespace std;

template<template<typename> typename R, typename T, typename Objective, typename Derivative>
requires(is_floating_point_v<T> &&
         is_invocable_r_v<T, Objective, R<T> const &> &&
         is_invocable_r_v<tuple<vector<T>, T>, Derivative, R<T> const &> &&
         is_base_of<DerivativeOptimizable<R, T>, R<T>>::value)
R<T> lineSearch(Objective && objective, Derivative && derivative, R<T> const & startingPoint,
                float stagnationTolerance, int maxLineSearchSteps, float maxStepLength,
                float tol, bool bfgs = false) {
    R<T> r = startingPoint.copy();
    T obj = invoke(forward<Objective>(objective), forward<R<T> const &>(startingPoint));
    cout << "LS starting objective " << obj << "\n";
    vector<T> deriv, deriv2, nextDeriv;
    T derivNorm, derivNorm2, nextDerivNorm;
    tie(deriv, derivNorm) = invoke(forward<Derivative>(derivative), forward<R<T> const &>(r));
    vector<T> stepDirection = deriv;
    T stepDirectionNorm = derivNorm;
    T maxStep = maxStepLength;
    int stepsSinceDerivUpdate = 0;
    T armijoLimit = 1E-4;
    T curvatureLimit = 0.9;
    VecUtil::Matrix<T> B = VecUtil::Matrix<T>::identity(deriv.size());
    while(true) {
        T minStep = 0;
        vector<T> steps;
        T bestStepLength = 0;
        bool satisfiedWolfeConditions = false;
        if(bfgs) {
            /*
             * eigenvalues, eigenvectors = scipy.linalg.eigh(B)
             * print("B eigs", eigenvalues)
             */
            T conditionNumberLimit = 10;
            if (false) { //numpy.max(numpy.abs(eigenvalues)) / numpy.min(numpy.abs(eigenvalues)) > conditionNumberLimit:
                cout << "condition number exceeded " << conditionNumberLimit << "restarting with B=I\n";
                B = VecUtil::Matrix<T>::identity(B.getNumRows());
            }
            stepDirection = B * deriv;
            stepDirectionNorm = VecUtil::norm(stepDirection);
            T angle = VecUtil::dot(stepDirection, deriv) / stepDirectionNorm / derivNorm;
            cout << "BFGS step angle with deriv " << angle << "\n";
            T angleLimit = 0.05;
            if(angle < angleLimit) {
                cout << "Angle limit exceeded, reseting B\n";
                B = VecUtil::Matrix<T>::identity(B.getNumRows());
                stepDirection = deriv;
                stepDirectionNorm = derivNorm;
            }
        } else {
            stepDirection = deriv;
            stepDirectionNorm = derivNorm;
        }
        for(int j = 0; j < maxLineSearchSteps; ++j) {
            vector<T> steps = Util::linspace(minStep, maxStep, 5);
            vector<T> objs;
            for(T h : steps) {
                R<T> newPoint = r.copy();
                newPoint.updateCoordinates(stepDirection, h);
                objs.push_back(invoke(forward<Objective>(objective), forward<R<T> const &>(newPoint)));
                if (h == 0) {
                    continue;
                }
                T t1 = -(objs.back() - obj) / h / (stepDirectionNorm * stepDirectionNorm);
                cout << "Armijo limit would be compared to " << t1 << " h: " << h << " obj: " << objs.back() << "\n";
                if (t1 >= armijoLimit) {
                    tie(deriv2, derivNorm2) = invoke(forward<Derivative>(derivative), forward<R<T> const &>(newPoint));
                    T t2 = VecUtil::dot(stepDirection, deriv2) / stepDirectionNorm;
                    cout << "***Curvature limit would be compared to " << t2 << " h: " << h << " obj: " << objs.back()
                         << "\n";
                    if (t2 < curvatureLimit) {
                        bestStepLength = h;
                        satisfiedWolfeConditions = true;
                        break;
                    }
                }
            }

            if(!satisfiedWolfeConditions) {
                int minIdx = Util::argmin(objs);
                if(abs(objs.back() - objs[0]) / abs(objs[0]) < 1E-2) {
                    break;
                }
                if(minIdx == 0) {
                    maxStep = steps[1];
                } else if(minIdx == steps.size() - 1) {
                    minStep = maxStep;
                    maxStep += 5 * (steps[steps.size()-1] - steps[steps.size()-2]);
                } else {
                    minStep = steps[minIdx - 1];
                    maxStep = steps[minIdx + 1];
                    bestStepLength = steps[minIdx];
                }
            } else {
                break;
            }
        }
        R<T> oldR = r.copy();
        r.updateCoordinates(stepDirection, bestStepLength);
        T maxStep = bestStepLength * 10;
        cout << "LS step size used "<<bestStepLength<<".  Wolfe conditions "<<(satisfiedWolfeConditions ? "satisfied":"not satisfied")<<"\n";
        T newObj = invoke(forward<Objective>(objective), forward<R<T> const &>(r));
        tie(nextDeriv, nextDerivNorm) = invoke(forward<Derivative>(derivative), forward<R<T> const &>(r));
        if(bfgs) {
            vector<T> y = VecUtil::subtract(nextDeriv, deriv);
            vector<T> s = VecUtil::subtract(r.getCoordinates(), oldR.getCoordinates());
            T rho = 1 / VecUtil::dot(y, s);
            cout << "rho " << rho << "\n";
            cout << "norm y " << VecUtil::norm(y) << "\n";
            cout << "norm s " << VecUtil::norm(s) << "\n";
            VecUtil::Matrix<T> updateMatrix(s.size(), s.size());
            updateMatrix.fill([rho, &y, &s](int i, int j) {
                if(i != j) {
                    return -rho*y[i]*s[j];
                } else {
                    return 1-rho*y[i]*s[j];
                }
            });
            VecUtil::Matrix<T> sMatrix(deriv.size());
            sMatrix.fill([rho, &s](int i, int j) {
                return rho*s[j]*s[i];
            });
            B = updateMatrix.transpose() * B * updateMatrix + sMatrix;
        }
        deriv = nextDeriv;
        derivNorm = nextDerivNorm;
        stepsSinceDerivUpdate += 1;
        if(stepsSinceDerivUpdate == 100) {
            stepsSinceDerivUpdate = 0;
        }
        cout << "LS step " << newObj << " " <<  derivNorm << "\n";
        if(newObj < 2) {
            bfgs = true;
        }
        if(abs((newObj - obj)/obj) < stagnationTolerance) {
            cout << "LS stagnated\n";
            break;
        }
        obj = newObj;
        if(derivNorm < tol) {
            cout << "LS converged\n";
            break;
        }
        if(bestStepLength == 0) {
            cout << "LS got stuck\n";
            break;
        }
    }
    return r;
}


#endif //CPP_LINESEARCH_H
