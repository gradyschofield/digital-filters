#include <iostream>
#include <vector>

#include<DiskRational.h>
#include<Util.h>

using namespace std;

template<typename T>
void lineSearch(vector<T> const & grid, vector<T> const & filter, int sampleRate, float stagnationTolerance,
                int maxLineSearchSteps, float maxStepLength, float tol, DiskRational<T> const & startingPoint,
                bool bfgs = false) {
    DiskRational<T> r = startingPoint.copy();
    T obj = Util::objective(grid, filter, sampleRate, r);
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
    vector<T> B(deriv.size()*deriv.size());
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
            if (false) { #numpy.max(numpy.abs(eigenvalues)) / numpy.min(numpy.abs(eigenvalues)) > conditionNumberLimit:
                print("condition number exceeded ", conditionNumberLimit, "restarting with B=I")
                B = numpy.identity(len(deriv))
            }
            stepDirection = matmul(B, deriv, deriv.size());
            stepDirectionNorm = VecUtil::norm(stepDirection);
            T angle = VecUtil::dot(stepDirection, deriv) / stepDirectionNorm / derivNorm;
            cout << "BFGS step angle with deriv" << angle << "\n";
            T angleLimit = 0.05;
            if(angle < angleLimit) {
                cout << "Angle limit exceeded, reseting B\n";
                VecUtil::fillIdentity(B);
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
                if(h == 0) {
                    continue;
                }
                T t1 = -(objs[-1] - obj) / h / (stepDirectionNorm * stepDirectionNorm);
                cout << "Armijo limit would be compared to " << t1 << "h: " <<  h << "obj: " <<  objs[-1] << "\n";
                if(t1 >= armijoLimit ) {
                    tie(deriv2, derivNorm2) = newPoint.derivative(grid, filter, sampleRate);
                    T t2 = VecUtil::dot(stepDirection, deriv2) / stepDirectionNorm;
                    cout<<"***Curvature limit would be compared to "<<t2<<"h: "<<h<<"obj: "<<objs.back()<<"\n";
                    if(t2 < curvatureLimit) {
                        bestStepLength = h;
                        satisfiedWolfeConditions = true;
                        break;
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
                } else if(minIdx == len(steps) - 1) {
                    minStep = maxStep;
                    maxStep += 5 * (steps[-1] - steps[-2]);
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
        cout << "LS step size used "<<bestStepLength<<".  Wolfe conditions "<<(satisfiedWolfeConditions?"satisfied":"not satisfied"<<"\n";
        T newDot = innerProduct(grid, filter, sampleRate, r);
        T newObj = Util::objective(grid, filter, sampleRate, r);
        tie(nextDeriv, nextDerivNorm) = r.derivative(grid, filter, sampleRate);
        if(bfgs) {
            vector<T> y = nextDeriv - deriv;
            vector<T> s = r.getCoordinates() - oldR.getCoordinates();
            T rho = 1 / y.dot(s);
            print("rho", rho)
            print("norm y", y.norm())
            print("norm s", s.norm())
            updateMatrix = [[(0 if i != j else 1) - rho*y[i]*s[j] for i in range(len(s))] for j in range(len(s))]
            sMatrix = [[rho*s[j]*s[i] for i in range(len(s))] for j in range(len(s))]
            B = numpy.matmul(
                    numpy.matmul(numpy.transpose(updateMatrix), B),
                    updateMatrix) + sMatrix
        }
        deriv = nextDeriv
        derivNorm = nextDerivNorm
        stepsSinceDerivUpdate += 1
        if stepsSinceDerivUpdate == 100:
        update(r, filter, sampleRate)
        stepsSinceDerivUpdate = 0
        print("LS step", newDot, newObj, derivNorm)
        if newObj < 2:
        bfgs = True
        runInfo.append(GDRunInfo(derivNorm, newDot, newObj))
        if abs((newObj - obj)/obj) < stagnationTolerance:
        print("LS stagnated")
        break
        obj = newObj
        if derivNorm < tol:
        print("LS converged")
        break
        if bestStepLength == 0:
        print("LS got stuck")
        break
    update(r, filter, sampleRate)
    x, y = r.plotData(0, 2000, 96000)
    plt.plot(x, y)
    plt.title("LS result zoomed")
    plt.show()
    x, y = r.plotData(0, 48000, 96000)
    plt.plot(x, y)
    plt.title("LS result")
    plt.show()
    return r, runInfo

        def geneticOptimizer(filter, maxGenerations, populationSize, cullSize,
        RationalType,
        numNumeratorReal, numNumeratorComplex,
        numDenominatorReal, numDenominatorComplex,
        startingPoint=None):
if startingPoint:
obj = Util.objective(grid, filter, sampleRate, startingPoint)
population = [(startingPoint, obj)]
else:
population = []
for i in range(populationSize - len(population)):
r = RationalType.create(numNumeratorReal, numNumeratorComplex, numDenominatorReal, numDenominatorComplex)
obj = Util.objective(grid, filter, sampleRate, r)
population.append((r, obj))

population.sort(key=lambda x: x[1])
del population[cullSize:-1]
print('best obj', population[0][1])
for j in range(maxGenerations):
for i in range(cullSize, populationSize):
r1 = population[random.randint(0, cullSize-1)][0]
r2 = population[random.randint(0, cullSize-1)][0]
r = r1.crossover(r2)
r.mutate(0.5, 0.01)
obj = Util.objective(grid, filter, sampleRate, r)
population.append((r, obj))
population.sort(key=lambda x: x[1])
del population[cullSize:-1]
print('best obj', population[0][1])
r = population[0][0]
x, y = r.plotData(0, 48000, 96000)
plt.plot(x, y)
plt.title("Best in generation")
plt.show()
return population[0][0]

int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
