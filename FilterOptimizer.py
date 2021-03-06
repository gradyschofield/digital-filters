import random

import numpy
import scipy.linalg
import math
import cmath
import functools
import random
import matplotlib.pyplot as plt

import Chebyshev
import DiskRational
import DiskRational2
import BSpline
import Rational
import VecUtil
import Util

class GDRunInfo:
    derivNorm = 0.0
    dot = 0.0
    objective = 0.0

    def __init__(self, derivNorm, dot, objective):
        self.derivNorm = derivNorm
        self.dot = dot
        self.objective = objective

def innerProduct(grid, filter, sampleRate, rational):
    stepIntegral = 0
    rationalIntegral = 0
    dotIntegral = 0
    for freq, filterValue in zip(grid, filter):
        x = cmath.exp(complex(0, -2 * math.pi * freq / sampleRate))
        y = abs(rational.evaluate(x))
        stepIntegral += filterValue * filterValue
        dotIntegral += y * filterValue
        rationalIntegral += y * y
    return dotIntegral / math.sqrt(stepIntegral * rationalIntegral)

def randomSearch(grid, filter, sampleRate, numSteps, numeratorReal, numeratorComplex, denominatorReal, denominatorComplex, RationalType):
    bestRational = None
    #dot = -2
    dot = 2**31
    for i in range(numSteps):
        r = RationalType.create(numeratorReal, numeratorComplex, denominatorReal, denominatorComplex)
        newDot = innerProduct(grid, filter, sampleRate, r)
        newDot = Util.objective(grid, filter, sampleRate, r)
        #if newDot > dot:
        if newDot < dot:
            bestRational = r
            dot = newDot
            print(i, dot)
    x, y = bestRational.plotData(0, 48000, 96000)
    plt.plot(x, y)
    plt.title("Random search result")
    plt.show()
    return bestRational

def newton(grid, filter, sampleRate, maxSteps, maxStepLength, tol, startingPoint, printEachStep = False, computeEigenvector = False):
    r = startingPoint.copy()
    runInfo = []
    obj = Util.objective(grid, filter, sampleRate, r)
    bestObj = obj
    print('newton starting objective', obj)
    deriv, derivNorm = r.derivative(grid, filter, sampleRate)
    for i in range(maxSteps):
        if derivNorm < tol:
            print("Newton converged")
            break
        deriv2 = r.secondDerivative(grid, filter, sampleRate)
        if computeEigenvector:
            eigenvalues, eigenvecors = scipy.linalg.eigh(deriv2)
            print("eigenvalues", eigenvalues)
            rc = list(r.denominatorCoef)
            rc.reverse()
            denominatorRoots = numpy.roots(rc)
            x = [z.real for z in denominatorRoots]
            y = [z.imag for z in denominatorRoots]
            plt.scatter(x, y, marker='o')
            rc = list(r.numeratorCoef)
            rc.reverse()
            numeratorRoots = numpy.roots(rc)
            x = [z.real for z in numeratorRoots]
            y = [z.imag for z in numeratorRoots]
            plt.scatter(x, y, marker='^')
            plt.title("LS partial result denominator roots")
            plt.show()
            x, y = r.plotData(0, 48000, 96000)
            plt.plot(x, y)
            plt.title("Newton result")
            plt.show()
        h = numpy.linalg.solve(deriv2, deriv)
        stepLength = VecUtil.norm(h)
        h = VecUtil.maxLength(h, maxStepLength)
        r.updateCoef(h, 1)
        newDot = innerProduct(grid, filter, sampleRate, r)
        obj = Util.objective(grid, filter, sampleRate, r)
        deriv, derivNorm = r.derivative(grid, filter, sampleRate)
        print("Newton step dot:", newDot, "obj:", obj, "derivNorm:", derivNorm, "step length:", stepLength)
        runInfo.append(GDRunInfo(derivNorm, newDot, obj))
        if bestObj > obj:
            bestObj = obj
            if printEachStep:
                x, y = r.plotData(0, 48000, 96000)
                plt.plot(x, y)
                plt.title("Newton result")
                plt.show()
    else:
        print('Newton reached maxSteps')
    x, y = r.plotData(0, 48000, 96000)
    plt.plot(x, y)
    plt.title("Newton result")
    plt.show()
    return r, runInfo

'''def gradientDescent(maxSteps, maxStepLength, tol, startingPoint):
    r = Rational(startingPoint.numeratorCoef, startingPoint.denominatorCoef)
    runInfo = []
    obj = r.objective(stepMin, stepMax, gain, filter, sampleRate)
    print('GD starting objective', obj)
    deriv, derivNorm = r.derivative(stepMin, stepMax, gain, filter, sampleRate)
    for i in range(maxSteps):
        maxLength(deriv, maxStepLength, rate)
        r.updateCoef(deriv, rate)
        newDot = r.innerProduct(stepMin, stepMax, gain, filter, sampleRate)
        obj = r.objective(stepMin, stepMax, gain, filter, sampleRate)
        deriv, derivNorm = r.derivative(stepMin, stepMax, gain, filter, sampleRate)
        print("GD step", newDot, obj, derivNorm)
        runInfo.append(GDRunInfo(derivNorm, newDot, obj))
        if derivNorm < tol:
            print("GD converged")
            break
    else:
        print("GD reached maxSteps")
    x, y = r.plotData(0, 48000, 96000)
    plt.plot(x, y)
    plt.title("GD result")
    plt.show()
    return r, runInfo'''

def update(r, sampleRate):
    '''deriv2 = r.secondDerivative(grid, filter, sampleRate)
    eigenvalues, eigenvecors = scipy.linalg.eigh(deriv2)
    print("eigenvalues", eigenvalues)'''
    x, y = r.plotData(0, sampleRate/2, sampleRate)
    plt.plot(x, y)
    plt.title("LS partial result")
    plt.show()
    numeratorRoots, denominatorRoots = r.getAllRoots()
    x = [z.real for z in denominatorRoots]
    y = [z.imag for z in denominatorRoots]
    plt.scatter(x, y, marker='o')
    x = [z.real for z in numeratorRoots]
    y = [z.imag for z in numeratorRoots]
    plt.scatter(x, y, marker='^')
    x = [math.cos(theta) for theta in numpy.linspace(0, 2 * math.pi, 100)]
    y = [math.sin(theta) for theta in numpy.linspace(0, 2 * math.pi, 100)]
    plt.plot(x, y)
    plt.title("LS partial result roots")
    plt.show()

    x = [z.real for z in denominatorRoots]
    y = [z.imag for z in denominatorRoots]
    plt.scatter(x, y, marker='o')
    plt.title("LS partial result denominator roots")
    plt.show()


def lineSearch(grid, filter, sampleRate, stagnationTolerance, maxLineSearchSteps, maxStepLength, tol, startingPoint, bfgs = False):
    r = startingPoint.copy()
    runInfo = []
    obj = Util.objective(grid, filter, sampleRate, r)
    print('LS starting objective', obj)
    deriv, derivNorm = r.derivative(grid, filter, sampleRate)
    stepDirection = deriv
    stepDirectionNorm = derivNorm
    maxStep = maxStepLength
    stepsSinceDerivUpdate = 0
    armijoLimit = 1E-4
    curvatureLimit = 0.9
    B = numpy.identity(len(deriv))
    while True:
        minStep = 0
        steps = []
        bestStepLength = 0
        satisfiedWolfeConditions = False
        if bfgs:
            eigenvalues, eigenvectors = scipy.linalg.eigh(B)
            print("B eigs", eigenvalues)
            conditionNumberLimit = 10
            if False:#numpy.max(numpy.abs(eigenvalues)) / numpy.min(numpy.abs(eigenvalues)) > conditionNumberLimit:
                print("condition number exceeded ", conditionNumberLimit, "restarting with B=I")
                B = numpy.identity(len(deriv))
            stepDirection = numpy.matmul(B, deriv)
            stepDirectionNorm = VecUtil.norm(stepDirection)
            angle = VecUtil.dot(stepDirection, deriv) / stepDirectionNorm / derivNorm
            print("BFGS step angle with deriv", angle)
            angleLimit = 0.05
            if angle < angleLimit:
                print("Angle limit exceeded, reseting B")
                B = numpy.identity(len(deriv))
                stepDirection = deriv
                stepDirectionNorm = derivNorm
        else:
            stepDirection = deriv
            stepDirectionNorm = derivNorm
        for j in range(maxLineSearchSteps):
            steps = numpy.linspace(minStep, maxStep, 5)
            objs = []
            for h in steps:
                newPoint = r.copy().updateCoef(stepDirection, h)
                objs.append(Util.objective(grid, filter, sampleRate, newPoint))
                if h == 0:
                    continue
                t1 = -(objs[-1] - obj) / h / stepDirectionNorm**2
                print("Armijo limit would be compared to ", t1, "h: ", h, "obj: ", objs[-1])
                if t1 >= armijoLimit:
                    deriv2, derivNorm2 = newPoint.derivative(grid, filter, sampleRate)
                    t2 = VecUtil.dot(stepDirection, deriv2) / stepDirectionNorm
                    print("***Curvature limit would be compared to ", t2, "h: ", h, "obj: ", objs[-1])
                    if t2 < curvatureLimit:
                        bestStepLength = h
                        satisfiedWolfeConditions = True
                        break

            if not satisfiedWolfeConditions:
                #objs = [objective(grid, filter, sampleRate, r.copy().updateCoef(deriv, rate)) for rate in steps]
                #print(objs)
                minIdx = numpy.argmin(objs)
                if abs(objs[-1] - objs[0])/abs(objs[0]) < 1E-2:
                    break
                if minIdx == 0:
                    maxStep = steps[1]
                elif minIdx == len(steps) - 1:
                    minStep = maxStep
                    maxStep += 5*(steps[-1] - steps[-2])
                else:
                    minStep = steps[minIdx-1]
                    maxStep = steps[minIdx+1]
                bestStepLength = steps[minIdx]
            else:
                break
        oldR = r.copy()
        r.updateCoef(stepDirection, bestStepLength)
        maxStep = bestStepLength * 10
        print('LS step size used', bestStepLength, 'Wolfe conditions', 'satisfied' if satisfiedWolfeConditions else 'not satisfied')
        newDot = innerProduct(grid, filter, sampleRate, r)
        newObj = Util.objective(grid, filter, sampleRate, r)
        nextDeriv, nextDerivNorm = r.derivative(grid, filter, sampleRate)
        if bfgs:
            y = [t1 - t2 for t1, t2 in zip(nextDeriv, deriv)]
            s = [t1 - t2 for t1, t2 in zip(r.getCoordinates(), oldR.getCoordinates())]
            rho = 1 / VecUtil.dot(y, s)
            print("rho", rho)
            print("norm y", VecUtil.norm(y))
            print("norm s", VecUtil.norm(s))
            updateMatrix = [[(0 if i != j else 1) - rho*y[i]*s[j] for i in range(len(s))] for j in range(len(s))]
            sMatrix = [[rho*s[j]*s[i] for i in range(len(s))] for j in range(len(s))]
            B = numpy.matmul(
                    numpy.matmul(numpy.transpose(updateMatrix), B),
                    updateMatrix) + sMatrix
        deriv = nextDeriv
        derivNorm = nextDerivNorm
        stepsSinceDerivUpdate += 1
        if stepsSinceDerivUpdate == 100:
            update(r, sampleRate)
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
    update(r, sampleRate)
    x, y = r.plotData(0, 2000, 96000)
    plt.plot(x, y)
    plt.title("LS result zoomed")
    plt.show()
    x, y = r.plotData(0, 48000, 96000)
    plt.plot(x, y)
    plt.title("LS result")
    plt.show()
    return r, runInfo

def geneticOptimizer(maxGenerations, populationSize, cullSize,
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


def testRoots(rational, log = True):
    if isinstance(rational, DiskRational.DiskRational):
        return True
    rc = list(rational.numeratorCoef)
    rc.reverse()
    numeratorRoots = numpy.roots(rc)
    rc = list(rational.denominatorCoef)
    rc.reverse()
    denominatorRoots = numpy.roots(rc)
    #b = functools.reduce(lambda b, z: abs(z) > 1 or b, numeratorRoots, False)
    b = functools.reduce(lambda b, z: abs(z) > 1 or b, denominatorRoots, False)
    if b:
        if log:
            print("Invalid roots")
            print('numerator roots:' ,abs(numeratorRoots))
            print('denominator roots:', abs(denominatorRoots))
        return False
    else:
        if log:
            print("Roots in unit disk")
        return True

if __name__ == '__main__':
    stepMin = 50
    stepMax = 100
    gain = 6
    gridPoints = 100
    isFilter = False
    sampleRate = 96000
    rate = 0.01
    degree = 500
    #rationalType = Rational.Rational
    rationalType = DiskRational.DiskRational
    numRealNumeratorRoots = 0
    numComplexNumeratorRoots = 20
    numRealDenominatorRoots = 0
    numComplexDenominatorRoots = 20

    grid = numpy.linspace(0, sampleRate/2, gridPoints)
    if False:
        c = Chebyshev.Chebyshev.projectStep(stepMin, stepMax, 0, sampleRate/2, gain, isFilter, degree)
        filterRoots = c.getAllRoots()
        x = [z.real for z in filterRoots]
        y = [z.imag for z in filterRoots]
        plt.scatter(x, y, marker='o')
        plt.title("Filter roots")
        plt.show()
        filter = c.evaluateOnScale(grid, 0, sampleRate/2)
    else:
        bspline = BSpline.BSpline(grid, BSpline.BSpline.getDefaultKnots())
        filter = bspline.getBasis(9)
    plt.plot(grid[1:100], filter[1:100])
    plt.title("Exact filter")
    plt.show()
    r = None
    while True:
        r = geneticOptimizer(5, 2000, 200, rationalType,
                         numRealNumeratorRoots, numComplexNumeratorRoots,
                         numRealDenominatorRoots, numComplexDenominatorRoots)
        '''r = randomSearch(grid, filter, sampleRate, 10, numRealNumeratorRoots,
                         numComplexNumeratorRoots, numRealDenominatorRoots,
                         numComplexDenominatorRoots, rationalType)'''
        #testRoots(r)
        #r, gdRunInfo = gradientDescent(20, 1E-7, r)
        r, lsRunInfo = lineSearch(grid, filter, sampleRate, 1E-5, 10, 10, 1E-7, r, False)
        #testRoots(r)
        #r, newtonRunInfo = newton(grid, filter, sampleRate, 100, 0, 1E-7, r, printEachStep=False, computeEigenvector=True)
        #testRoots(r)
        x, y = r.plotData(0, 48000, 96000, gridPoints)
        residual = [t1 - t2 for t1, t2 in zip(y, filter)]
        plt.plot(grid, residual)
        plt.title('Residual')
        plt.show()
