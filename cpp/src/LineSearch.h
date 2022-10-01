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

template<template<typename> typename Optimizable, typename FieldElement, typename ObjectiveFunction, typename DerivativeFunction>
requires(is_floating_point_v<FieldElement> &&
         is_invocable_r_v<FieldElement, ObjectiveFunction, Optimizable<FieldElement> const &> &&
         is_invocable_r_v<tuple<vector<FieldElement>, FieldElement>, DerivativeFunction, Optimizable<FieldElement> const &> &&
         is_base_of<DerivativeOptimizable<Optimizable, FieldElement>, Optimizable<FieldElement>>::value)
Optimizable<FieldElement> lineSearch(ObjectiveFunction && computeObjective, DerivativeFunction && computeDerivative, Optimizable<FieldElement> const & startingPoint,
                                     float stagnationTolerance, int maxLineSearchSteps, float maxStepLength,
                                     float tol, bool bfgs = false) {
    Optimizable<FieldElement> r = startingPoint.copy();
    FieldElement obj = invoke(forward<ObjectiveFunction>(computeObjective), forward<Optimizable<FieldElement> const &>(startingPoint));
    cout << "LS starting objective " << obj << "\n";
    vector<FieldElement> deriv, deriv2, nextDeriv;
    FieldElement derivNorm, derivNorm2, nextDerivNorm;
    tie(deriv, derivNorm) = invoke(forward<DerivativeFunction>(computeDerivative), forward<Optimizable<FieldElement> const &>(r));
    vector<FieldElement> stepDirection = deriv;
    FieldElement stepDirectionNorm = derivNorm;
    FieldElement maxStep = maxStepLength;
    int stepsSinceDerivUpdate = 0;
    FieldElement armijoLimit = 1E-4;
    FieldElement curvatureLimit = 0.9;
    VecUtil::Matrix<FieldElement> B = VecUtil::Matrix<FieldElement>::identity(deriv.size());
    while(true) {
        FieldElement minStep = 0;
        vector<FieldElement> steps;
        FieldElement bestStepLength = 0;
        bool satisfiedWolfeConditions = false;
        if(bfgs) {
            /*
             * eigenvalues, eigenvectors = scipy.linalg.eigh(B)
             * print("B eigs", eigenvalues)
             */
            FieldElement conditionNumberLimit = 10;
            if (false) { //numpy.max(numpy.abs(eigenvalues)) / numpy.min(numpy.abs(eigenvalues)) > conditionNumberLimit:
                cout << "condition number exceeded " << conditionNumberLimit << "restarting with B=I\n";
                B = VecUtil::Matrix<FieldElement>::identity(B.getNumRows());
            }
            stepDirection = B * deriv;
            stepDirectionNorm = VecUtil::norm(stepDirection);
            FieldElement angle = VecUtil::dot(stepDirection, deriv) / stepDirectionNorm / derivNorm;
            cout << "BFGS step angle with deriv " << angle << "\n";
            FieldElement angleLimit = 0.05;
            if(angle < angleLimit) {
                cout << "Angle limit exceeded, reseting B\n";
                B = VecUtil::Matrix<FieldElement>::identity(B.getNumRows());
                stepDirection = deriv;
                stepDirectionNorm = derivNorm;
            }
        } else {
            stepDirection = deriv;
            stepDirectionNorm = derivNorm;
        }
        for(int j = 0; j < maxLineSearchSteps; ++j) {
            vector<FieldElement> steps = Util::linspace(minStep, maxStep, 5);
            vector<FieldElement> objs;
            for(FieldElement h : steps) {
                Optimizable<FieldElement> newPoint = r.copy();
                newPoint.updateCoordinates(stepDirection, h);
                objs.push_back(invoke(forward<ObjectiveFunction>(computeObjective), forward<Optimizable<FieldElement> const &>(newPoint)));
                if (h == 0) {
                    continue;
                }
                FieldElement t1 = -(objs.back() - obj) / h / (stepDirectionNorm * stepDirectionNorm);
                cout << "Armijo limit would be compared to " << t1 << " h: " << h << " obj: " << objs.back() << "\n";
                if (t1 >= armijoLimit) {
                    tie(deriv2, derivNorm2) = invoke(forward<DerivativeFunction>(computeDerivative), forward<Optimizable<FieldElement> const &>(newPoint));
                    FieldElement t2 = VecUtil::dot(stepDirection, deriv2) / stepDirectionNorm;
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
        Optimizable<FieldElement> oldR = r.copy();
        r.updateCoordinates(stepDirection, bestStepLength);
        FieldElement maxStep = bestStepLength * 10;
        cout << "LS step size used "<<bestStepLength<<".  Wolfe conditions "<<(satisfiedWolfeConditions ? "satisfied":"not satisfied")<<"\n";
        FieldElement newObj = invoke(forward<ObjectiveFunction>(computeObjective), forward<Optimizable<FieldElement> const &>(r));
        tie(nextDeriv, nextDerivNorm) = invoke(forward<DerivativeFunction>(computeDerivative), forward<Optimizable<FieldElement> const &>(r));
        if(bfgs) {
            vector<FieldElement> y = VecUtil::subtract(nextDeriv, deriv);
            vector<FieldElement> s = VecUtil::subtract(r.getCoordinates(), oldR.getCoordinates());
            FieldElement rho = 1 / VecUtil::dot(y, s);
            cout << "rho " << rho << "\n";
            cout << "norm y " << VecUtil::norm(y) << "\n";
            cout << "norm s " << VecUtil::norm(s) << "\n";
            VecUtil::Matrix<FieldElement> updateMatrix(s.size(), s.size());
            updateMatrix.fill([rho, &y, &s](int i, int j) {
                if(i != j) {
                    return -rho*y[i]*s[j];
                } else {
                    return 1-rho*y[i]*s[j];
                }
            });
            VecUtil::Matrix<FieldElement> sMatrix(deriv.size());
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
