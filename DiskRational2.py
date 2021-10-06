import cmath
import math
import random
import functools
import numpy

import Util
import BaseRational
import BSpline

class DiskRational2(BaseRational.BaseRational):
    numeratorInnerRadius = 0.5
    numeratorOuterRadius = 0.99999
    denominatorInnerRadius = 0
    denominatorOuterRadius = 0.1

    numeratorRoots = []
    denominatorRoots = []

    numeratorX = []
    numeratorThetaX = []

    denominatorX = []
    denominatorThetaX = []

    def __init__(self, numeratorRoots, denominatorRoots):
        self.numeratorRoots = list(numeratorRoots)
        self.denominatorRoots = list(denominatorRoots)
        self.computeParamsFromRoots()

    @staticmethod
    def computeParamsFromRootsHelper(numeratorRoots, denominatorRoots):
        numeratorX = []
        numeratorThetaX = []
        for r in numeratorRoots:
            length = abs(r)
            t = (length - DiskRational2.numeratorInnerRadius)/(DiskRational2.numeratorOuterRadius - DiskRational2.numeratorInnerRadius)
            numeratorX.append(-math.log(1/t - 1))
            theta = math.atan2(abs(r.imag), r.real)
            numeratorThetaX.append(-math.log(math.pi/theta - 1))
        denominatorX = []
        denominatorThetaX = []
        for r in denominatorRoots:
            length = abs(r)
            t = (length - DiskRational2.denominatorInnerRadius)/(DiskRational2.denominatorOuterRadius - DiskRational2.denominatorInnerRadius)
            denominatorX.append(-math.log(1/t - 1))
            theta = math.atan2(abs(r.imag), r.real)
            denominatorThetaX.append(-math.log(math.pi/theta - 1))
        return numeratorX, numeratorThetaX, denominatorX, denominatorThetaX

    def computeParamsFromRoots(self):
        self.numeratorX, \
        self.numeratorThetaX, \
        self.denominatorX, \
        self.denominatorThetaX = self.computeParamsFromRootsHelper(self.numeratorRoots, self.denominatorRoots)

    def getRoots(self):
        return list(self.numeratorRoots), list(self.denominatorRoots)

    def getAllRoots(self):
        numeratorRoots = list(self.numeratorRoots)
        numeratorRoots.extend([z.conjugate() for z in self.numeratorRoots])
        denominatorRoots = list(self.denominatorRoots)
        denominatorRoots.extend([z.conjugate() for z in self.denominatorRoots])
        return numeratorRoots, denominatorRoots

    @staticmethod
    def computeRootsFromParamsHelper(x, thetaX, innerRadius, outerRadius):
        r = []
        for x, thetaX in zip(x, thetaX):
            theta = math.pi / (1+math.exp(-thetaX))
            z = (innerRadius + (outerRadius-innerRadius)/(1+math.exp(-x))) * cmath.exp(complex(0,theta))
            r.append(z)
        return r

    def computeRootsFromParams(self):
        self.numeratorRoots = self.computeRootsFromParamsHelper( self.numeratorX, self.numeratorThetaX,
                                                                 DiskRational2.numeratorInnerRadius, DiskRational2.numeratorOuterRadius)
        self.denominatorRoots = self.computeRootsFromParamsHelper( self.denominatorX, self.denominatorThetaX,
                                                                   DiskRational2.denominatorInnerRadius, DiskRational2.denominatorOuterRadius)

    @staticmethod
    def create(dummy1, numNumeratorRoots, dummy2, numDenominatorRoots):
        numeratorRoots = Util.filterRootsHalfAnnulus(numNumeratorRoots, DiskRational2.numeratorInnerRadius, DiskRational2.numeratorOuterRadius)
        denominatorRoots = Util.filterRootsHalfAnnulus(numDenominatorRoots, DiskRational2.denominatorInnerRadius, DiskRational2.denominatorOuterRadius)
        return DiskRational2(numeratorRoots, denominatorRoots)

    def copy(self):
        return DiskRational2(self.numeratorRoots, self.denominatorRoots)

    def evaluate(self, x):
        n = 1
        d = 1
        for r in self.numeratorRoots:
            n *= (x - r)
            n *= (x - r.conjugate())
        for r in self.denominatorRoots:
            d *= (x - r)
            d *= (x - r.conjugate())
        return n / d

    def evaluateDenominator(self, x):
        d = 1
        for r in self.denominatorRoots:
            d *= (x - r)
            d *= (x - r.conjugate())
        return d

    def updateCoef(self, deriv, rate):
        k = 0
        for i in range(len(self.numeratorX)):
            self.numeratorX[i] -= rate * deriv[k]
            self.numeratorThetaX[i] -= rate * deriv[k+1]
            k += 2
        for i in range(len(self.denominatorX)):
            self.denominatorX[i] -= rate * deriv[k]
            self.denominatorThetaX[i] -= rate * deriv[k+1]
            k += 2
        self.computeRootsFromParams()
        return self

    def getCoordinates(self):
        ret = []
        ret.extend(self.numeratorX)
        ret.extend(self.numeratorThetaX)
        ret.extend(self.denominatorX)
        ret.extend(self.denominatorThetaX)
        return ret

    def derivative(self, grid, filter, sampleRate):
        deriv = [0 for i in range(len(self.numeratorX) +
                                  len(self.numeratorThetaX) +
                                  len(self.denominatorX) +
                                  len(self.denominatorThetaX))]
        for freq, filterValue in zip(grid, filter):
            z = cmath.exp(complex(0, -2 * math.pi * freq / sampleRate))
            y = self.evaluate(z)
            d = self.evaluateDenominator(z)
            f = abs(y)
            prefactor = 2 * (f - filterValue) / f
            k = 0
            for x, thetaX in zip(self.numeratorX, self.numeratorThetaX):
                theta = math.pi / (1 + math.exp(-thetaX))
                annulusWidth = DiskRational2.numeratorOuterRadius - DiskRational2.numeratorInnerRadius
                radius = DiskRational2.numeratorInnerRadius + annulusWidth / (1 + math.exp(-x))
                r = cmath.exp(complex(0, theta)) * radius
                z2 = z - r
                z2rconj = z - r.conjugate()
                dX = annulusWidth * math.exp(-x)/(1+math.exp(-x))**2 * cmath.exp(complex(0, theta))
                dTheta = radius * complex(0, math.pi) * cmath.exp(complex(0, theta)) * math.exp(-thetaX)/(1+math.exp(-thetaX))**2
                deriv[k] += prefactor * (y.conjugate() * y * (-dX / z2 - dX.conjugate() / z2rconj)).real
                deriv[k+1] += prefactor * (y.conjugate() * y * (-dTheta / z2 - dTheta.conjugate() / z2rconj)).real
                k += 2
            for x, thetaX in zip(self.denominatorX, self.denominatorThetaX):
                theta = math.pi / (1 + math.exp(-thetaX))
                annulusWidth = DiskRational2.denominatorOuterRadius - DiskRational2.denominatorInnerRadius
                radius = DiskRational2.denominatorInnerRadius + annulusWidth / (1 + math.exp(-x))
                r = cmath.exp(complex(0, theta)) * radius
                z2 = z - r
                z2rconj = z - r.conjugate()
                dX = annulusWidth * math.exp(-x) / (1 + math.exp(-x)) ** 2 * cmath.exp(complex(0, theta))
                dTheta = radius * complex(0, math.pi) * cmath.exp(complex(0, theta)) * math.exp(-thetaX) / (1 + math.exp(-thetaX)) ** 2
                deriv[k] += prefactor * (y.conjugate() * y * (dX / z2 + dX.conjugate() / z2rconj)).real
                deriv[k + 1] += prefactor * (y.conjugate() * y * (dTheta / z2 + dTheta.conjugate() / z2rconj)).real
                k += 2
        numGridPoints = len(grid)
        deriv = [t / numGridPoints for t in deriv]
        norm = functools.reduce(lambda n, t: n+t*t, deriv, 0)
        return deriv, math.sqrt(norm)

    def mutate(self, chanceOfMutation, mutationSize):
        newCoef = []
        for x in self.numeratorX:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
            newCoef.append(x)
        self.numeratorX = newCoef
        newCoef = []
        for x in self.numeratorThetaX:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
            newCoef.append(x)
        self.numeratorThetaX = newCoef
        newCoef = []
        for x in self.denominatorX:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
            newCoef.append(x)
        self.denominatorX = newCoef
        newCoef = []
        for x in self.denominatorThetaX:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
            newCoef.append(x)
        self.denominatorThetaX = newCoef
        self.computeRootsFromParams()


    def crossover(self, other):
        newNumeratorRoots = []
        for x, y in zip(self.numeratorRoots, other.numeratorRoots):
            newNumeratorRoots.append(x if random.uniform(0, 1) < 0.5 else y)
        newDenominatorRoots = []
        for x, y in zip(self.denominatorRoots, other.denominatorRoots):
            newDenominatorRoots.append(x if random.uniform(0, 1) < 0.5 else y)
        return DiskRational2(newNumeratorRoots, newDenominatorRoots)

if __name__ == '__main__':
    gridPoints = 200
    sampleRate = 96000

    grid = numpy.linspace(0, sampleRate/2, gridPoints)
    bspline = BSpline.BSpline(grid, BSpline.BSpline.getDefaultKnots())
    filter = bspline.getBasis(8)

    r = DiskRational2.create(0, 2, 0, 2)
    deriv, derivNorm = r.derivative(grid, filter, sampleRate)
    h = 0.001
    r2 = r.copy()
    r2.denominatorThetaX[0] += h
    r2.computeRootsFromParams()
    d = (Util.objective(grid, filter, sampleRate, r2) - Util.objective(grid, filter, sampleRate, r)) / h
    print(len(r.numeratorX))
    print(len(r.numeratorThetaX))
    print(len(r.denominatorX))
    print(len(r.denominatorThetaX))
    print(d)
    print(deriv[4])
    print(deriv)
