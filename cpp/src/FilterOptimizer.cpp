#include <iostream>
#include <fstream>
#include <vector>

#include<DiskRational.h>
#include<Util.h>
#include<VecUtil.h>
#include<BSpline.h>

using namespace std;


/*
void update(r, filter, sampleRate) {
    '''deriv2 = r.secondDerivative(grid, filter, sampleRate)
    eigenvalues, eigenvecors = scipy.linalg.eigh(deriv2)
    print("eigenvalues", eigenvalues)'''
    x, y = r.plotData(0, sampleRate/2, sampleRate)
    plt.plot(x, y)
    plt.title("LS partial result")
    plt.show()
    x, y = r.plotData(0, 48000, 96000, len(filter))
    residual = [t1 - t2 for t1, t2 in zip(y, filter)]
    plt.plot(grid, residual)
    plt.plot(grid, [0 for x in grid])
    plt.title('LS partial result residual')
    plt.show()
    plotRoots(r, "LS partial result roots")
    plotDenominatorRoots(r, "LS partial result denominator roots")
}
 */


template<typename T>
DiskRational<T> lineSearch(vector<T> const & grid, vector<T> const & filter, int sampleRate, float stagnationTolerance,
                           int maxLineSearchSteps, float maxStepLength, float tol, DiskRational<T> const & startingPoint,
                           bool bfgs = false) {
    DiskRational<T> r = startingPoint.copy();
    T obj = Util::objective(grid, filter, sampleRate, startingPoint);
    cout << "LS starting objective " << obj << "\n";
    vector<T> deriv, deriv2, nextDeriv;
    T derivNorm, derivNorm2, nextDerivNorm;
    tie(deriv, derivNorm) = r.derivative(grid, filter, sampleRate);
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
                DiskRational<T> newPoint = r.copy();
                newPoint.updateCoef(stepDirection, h);
                objs.push_back(Util::objective(grid, filter, sampleRate, newPoint));
                if (h == 0) {
                    continue;
                }
                T t1 = -(objs.back() - obj) / h / (stepDirectionNorm * stepDirectionNorm);
                cout << "Armijo limit would be compared to " << t1 << " h: " << h << " obj: " << objs.back() << "\n";
                if (t1 >= armijoLimit) {
                    tie(deriv2, derivNorm2) = newPoint.derivative(grid, filter, sampleRate);
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
                //objs =[objective(grid, filter, sampleRate, r.copy().updateCoef(deriv, rate)) for rate in steps]
                //print(objs)
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
        DiskRational<T> oldR = r.copy();
        r.updateCoef(stepDirection, bestStepLength);
        T maxStep = bestStepLength * 10;
        cout << "LS step size used "<<bestStepLength<<".  Wolfe conditions "<<(satisfiedWolfeConditions ? "satisfied":"not satisfied")<<"\n";
        T newDot = Util::innerProduct(grid, filter, sampleRate, r);
        T newObj = Util::objective(grid, filter, sampleRate, r);
        tie(nextDeriv, nextDerivNorm) = r.derivative(grid, filter, sampleRate);
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
            //update(r, filter, sampleRate);
            stepsSinceDerivUpdate = 0;
        }
        cout << "LS step " << newDot << " " <<  newObj << " " <<  derivNorm << "\n";
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
    //update(r, filter, sampleRate);
    /*
    x, y = r.plotData(0, 2000, 96000)
    plt.plot(x, y)
    plt.title("LS result zoomed")
    plt.show()
    x, y = r.plotData(0, 48000, 96000)
    plt.plot(x, y)
    plt.title("LS result")
    plt.show()
     */
    return r;
}

template<typename T>
DiskRational<T> geneticOptimizer(vector<T> const & grid, vector<T> const & filter,
                                 int sampleRate, int maxGenerations, int populationSize,
                                 int cullSize, int numNumeratorRoots, int numDenominatorRoots) {
    vector<pair<DiskRational<T>, T>> population;
    for(int i = 0; i < populationSize - population.size(); ++i) {
        DiskRational<T> r = DiskRational<T>::create(numNumeratorRoots, numDenominatorRoots);
        T obj = Util::objective(grid, filter, sampleRate, r);
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
            DiskRational<T> const & r1 = population[uniformIntDistribution(Util::rand::randomDevice)].first;
            DiskRational<T> const & r2 = population[uniformIntDistribution(Util::rand::randomDevice)].first;
            DiskRational<T> r = r1.crossover(r2);
            r.mutate(0.5, 0.01);
            T obj = Util::objective(grid, filter, sampleRate, r);
            population.emplace_back(move(r), obj);
        }
        sort(begin(population), end(population), [](auto & x, auto & y){
            return x.second < y.second;
        });
        population.erase(begin(population) + cullSize, end(population));
        cout << "best obj " << population[0].second << "\n";
        /*
        r = population[0][0]
        x, y = r.plotData(0, 48000, 96000)
        plt.plot(x, y)
        plt.title("Best in generation")
        plt.show()
         */
    }
    T obj = Util::objective(grid, filter, sampleRate, population[0].first);
    cout << "obj obj: " << obj << "\n";
    return move(population[0].first);
}

int main() {

    typedef float T;

    int sampleRate = 96000;
    int numGridPoints = 200;
    int numNumeratorRoots = 20;
    int numDenominatorRoots = 20;

    vector<T> grid = Util::linspace<T>(0, sampleRate/2, numGridPoints);
    BSpline<T> splineBasis(grid, BSpline<T>::getDefaultKnots());

    vector<T> filterFunc = splineBasis.getBasis(7);

    DiskRational<T> rational = geneticOptimizer(grid, filterFunc, sampleRate, 5, 2000, 200, numNumeratorRoots, numDenominatorRoots);
    ofstream rootTalk("rootTalk");
    while(true) {
        bool insertNumerator = true;
        rational = lineSearch(grid, filterFunc, sampleRate, 1E-5, 10, 10, 1E-7, rational, true);
        T minObj = numeric_limits<T>::max();
        complex<T> minRoot;
        for (T theta: Util::linspace<T>(0, M_PI, 20)) {
            for (T radius: Util::linspace<T>(0.01, insertNumerator ? 1.1 : 0.5, 20)) {
                complex<T> root = Util::toComplex(radius, theta);
                T obj;
                if (insertNumerator) {
                    obj = Util::objective(grid, filterFunc, sampleRate, rational, {root});
                } else {
                    obj = Util::objective(grid, filterFunc, sampleRate, rational, {}, {root});
                }
                if (obj < minObj) {
                    minObj = obj;
                    minRoot = root;
                }
            }
        }
        if (insertNumerator) {
            rational.incorporateRoots({minRoot}, {});
        } else {
            rational.incorporateRoots({}, {minRoot});
        }
        insertNumerator = !insertNumerator;
        rootTalk << minObj << endl;
    }

    return 0;
}
