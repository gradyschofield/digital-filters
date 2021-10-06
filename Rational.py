import math
import cmath
import matplotlib.pyplot as plt
import numpy
import random

import Util
import BaseRational

class Rational(BaseRational.BaseRational):
    numeratorCoef = []
    denominatorCoef = []

    def __init__(self, numeratorCoef, denominatorCoef):
        self.numeratorCoef = list(numeratorCoef)
        self.denominatorCoef = list(denominatorCoef)

    @staticmethod
    def create(numRealNumeratorRoots, numComplexNumeratorRoots, numRealDenominatorRoots, numComplexDenominatorRoots):
        numerator = Util.filterPolynomialFromRoots(numRealNumeratorRoots, numComplexNumeratorRoots, 1.1)
        denominator = Util.filterPolynomialFromRoots(numRealDenominatorRoots, numComplexDenominatorRoots, 0.99999)
        return Rational(numerator.coefficients, denominator.coefficients)

    @staticmethod
    def createWithCircularRootPattern(numRoots, numeratorRadius, denominatorRadius):
        numerator = Util.filterPolynomialWithCircularRootPattern(numRoots, numeratorRadius, 0, 1.0)
        denominator = Util.filterPolynomialWithCircularRootPattern(numRoots, denominatorRadius, 0.0, 1)
        return Rational(numerator.coefficients, denominator.coefficients)
        pass

    def copy(self):
        return Rational(self.numeratorCoef, self.denominatorCoef)

    def evaluate(self, x):
        n = 0
        d = 0
        for i in range(len(self.numeratorCoef)):
            n += x**i * self.numeratorCoef[i]
        for i in range(len(self.denominatorCoef)):
            d += x**i * self.denominatorCoef[i]
        return n / d

    def evaluateDenominator(self, x):
        d = 0
        for i in range(len(self.denominatorCoef)):
            d += x**i * self.denominatorCoef[i]
        return d

    def evaluatePartialNumerators(self, x):
        partialNumerators = []
        n = 0
        for i in range(len(self.numeratorCoef)):
            tmp = x**i * self.numeratorCoef[i]
            partialNumerators.append(-tmp)
            n += tmp
        return [x + n for x in partialNumerators]

    def getRoots(self):
        rc = list(self.numeratorCoef)
        rc.reverse()
        numeratorRoots = numpy.roots(rc)
        rc = list(self.denominatorCoef)
        rc.reverse()
        denominatorRoots = numpy.roots(rc)
        return numeratorRoots, denominatorRoots

    def dr(self, x):
        Dr = []
        Re_rDr = []
        for i in self.numeratorCoef:
            Dr.append(0.0)
            Re_rDr.append(0.0)
        for i in self.denominatorCoef:
            Dr.append(0.0)
            Re_rDr.append(0.0)
        y = self.evaluate(x)
        d = self.evaluateDenominator(x)
        k = 0
        for i in range(len(self.numeratorCoef)):
            Re_rDr[k] += (y.conjugate() * x**i / d).real
            Dr[k] += x**i / d
            k += 1
        for i in range(len(self.denominatorCoef)):
            Re_rDr[k] += -(y.conjugate() * y * x**i / d).real
            Dr[k] += -y * x**i / d
            k += 1
        return Re_rDr, Dr

    def derivative(self, grid, filter, sampleRate):
        deriv = []
        for x in self.numeratorCoef:
            deriv.append(0.0)
        for x in self.denominatorCoef:
            deriv.append(0.0)
        for freq, filterValue in zip(grid, filter):
            x = cmath.exp(complex(0, -2 * math.pi * freq / sampleRate))
            y = self.evaluate(x)
            f = abs(y)
            prefactor = 2 * (f - filterValue) / f
            d = self.evaluateDenominator(x)
            k = 0
            for i in range(len(self.numeratorCoef)):
                deriv[k] += prefactor * (y.conjugate() * x ** i / d).real
                k += 1
            for i in range(len(self.denominatorCoef)):
                deriv[k] += -prefactor * (y.conjugate() * y * x ** i / d).real
                k += 1
        numGridPoints = len(grid)
        deriv = [x / numGridPoints for x in deriv]
        norm = 0
        for x in deriv:
            norm += x ** 2
        return deriv, norm

    def getCoordinates(self):
        ret = []
        ret.extend(self.numeratorCoef)
        ret.extend(self.denominatorCoef)
        return ret

    def secondDerivative(self, grid, filter, sampleRate):
        size = len(self.numeratorCoef) + len(self.denominatorCoef)
        deriv2 = [[0 for j in range(size)] for i in range(size)]
        for freq, filterValue in zip(grid, filter):
            x = cmath.exp(complex(0, -2 * math.pi * freq / sampleRate))
            y = self.evaluate(x)
            f = abs(y)
            Re_rDr, Dr = self.dr(x)
            prefactor1 = 0.5 * filterValue / (y.conjugate() * y).real ** 1.5
            for i in range(len(Re_rDr)):
                for j in range(len(Re_rDr)):
                    deriv2[i][j] += prefactor1 * 2 * Re_rDr[i] * 2 * Re_rDr[j]
            prefactor2 = 1 - filterValue / f
            for i in range(len(Dr)):
                for j in range(len(Dr)):
                    deriv2[i][j] += prefactor2 * 2 * (Dr[i].conjugate() * Dr[j]).real
            monomialPowers = []
            for i in range(len(self.numeratorCoef)):
                monomialPowers.append(i)
            for i in range(len(self.denominatorCoef)):
                monomialPowers.append(i)
            d = self.evaluateDenominator(x)
            for i in range(len(Dr)):
                for j in range(len(Dr)):
                    ip = monomialPowers[i]
                    jp = monomialPowers[j]
                    if i < len(self.numeratorCoef) and j < len(self.numeratorCoef):
                        pass
                    elif i >= len(self.numeratorCoef) and j >= len(self.numeratorCoef):
                        deriv2[i][j] += prefactor2 * 2 * (y.conjugate() * y * x ** (ip + jp) / d ** 2).real
                    else:
                        deriv2[i][j] += -prefactor2 * 2 * (y.conjugate() * x ** (ip + jp) / d ** 2).real
        numGridPoints = len(grid)
        deriv2 = [[x / numGridPoints for x in row] for row in deriv2]
        '''norm = 0
        for i in range(len(Dr)):
            for j in range(len(Dr)):
                if i == j:
                    pass
                else:
                    norm += (deriv2[i][j] - deriv2[j][i])**2
        print("symmetry metric: ", math.sqrt(norm))'''
        '''norm = 0
        for i in range(len(Dr)):
            for j in range(len(Dr)):
                if i == j:
                    norm += (1-deriv2[i][j])**2
                else:
                    norm += deriv2[i][j]**2
        print("nearness to identity: ", math.sqrt(norm))'''
        return deriv2

    def updateCoef(self, deriv, rate):
        k = 0
        for i in range(len(self.numeratorCoef)):
            self.numeratorCoef[i] -= rate * deriv[k]
            k += 1
        for i in range(len(self.denominatorCoef)):
            self.denominatorCoef[i] -= rate * deriv[k]
            k += 1
        return self

    def mutate(self, chanceOfMutation, mutationSize):
        newCoef = []
        for x in self.numeratorCoef:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1+random.uniform(0, mutationSize)
            newCoef.append(x)
        self.numeratorCoef = newCoef
        newCoef = []
        for x in self.denominatorCoef:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1+random.uniform(0, mutationSize)
            newCoef.append(x)
        self.denominatorCoef = newCoef

    def crossover(self, other):
        newNumeratorCoef = []
        for x, y in zip(self.numeratorCoef, other.numeratorCoef):
            newNumeratorCoef.append(x if random.uniform(0,1) < 0.5 else y)
        newDenominatorCoef = []
        for x, y in zip(self.denominatorCoef, other.denominatorCoef):
            newDenominatorCoef.append(x if random.uniform(0,1) < 0.5 else y)
        return Rational(newNumeratorCoef, newDenominatorCoef)


if __name__ == '__main__':
    sampleRate = 96000
    gridPoints = 1000
    r = Rational.createWithCircularRootPattern(20, 0.8, 0.81)
    grid = numpy.linspace(0, sampleRate/2, gridPoints)

    x, y = r.plotData(0, 48000, 96000)
    plt.plot(x, y)
    plt.title("Circ")
    plt.show()
