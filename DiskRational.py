import cmath
import math
import random
import functools
import numpy

import Util
import BaseRational
import BSpline

class DiskRational(BaseRational.BaseRational):
    maxNumeratorRadius = 1.2
    maxDenominatorRadius = 0.2
    numeratorCoef = []

    realDenominatorX = []
    complexDenominatorX = []
    complexDenominatorTheta = []

    def __init__(self, numeratorCoef, realDenominatorRoots, complexDenominatorRoots):
        self.numeratorCoef = list(numeratorCoef)
        self.realDenominatorRoots = list(realDenominatorRoots)
        self.complexDenominatorRoots = list(complexDenominatorRoots)
        self.computeParamsFromRoots()
        self.fixedNumeratorRoots = []
        self.fixedDenominatorRoots = []

    @staticmethod
    def computeParamsFromRootsHelper(rr, cr, maxRadius):
        realX = []
        for r in rr:
            realX.append(-math.log(2/(r/maxRadius+1) - 1))
        complexX = []
        complexTheta = []
        for r in cr:
            complexX.append(-math.log(maxRadius/abs(r) - 1))
            complexTheta.append(math.atan2(r.imag, r.real))
        return realX, complexX, complexTheta

    def computeParamsFromRoots(self):
        self.realDenominatorX, \
        self.complexDenominatorX, \
        self.complexDenominatorTheta = self.computeParamsFromRootsHelper(self.realDenominatorRoots, self.complexDenominatorRoots, DiskRational.maxDenominatorRadius)

    def getRoots(self):
        rc = list(self.numeratorCoef)
        rc.reverse()
        numeratorRoots = numpy.roots(rc)
        denominatorRoots = [complex(x, 0) for x in self.realDenominatorRoots]
        denominatorRoots.extend(self.complexDenominatorRoots)
        return numeratorRoots, denominatorRoots

    def getAllRoots(self):
        rc = list(self.numeratorCoef)
        rc.reverse()
        numeratorRoots = numpy.roots(rc)
        denominatorRoots = [complex(x, 0) for x in self.realDenominatorRoots]
        denominatorRoots.extend(self.complexDenominatorRoots)
        denominatorRoots.extend([z.conjugate() for z in self.complexDenominatorRoots])
        return numeratorRoots, denominatorRoots

    @staticmethod
    def computeRootsFromParamsHelper(rx, cx, ctheta, maxRadius):
        rr = []
        for x in rx:
            rr.append(maxRadius * (1-math.exp(-x)/(1+math.exp(-x))))
        cr = []
        for x, theta in zip(cx, ctheta):
            z = maxRadius/(1+math.exp(-x)) * cmath.exp(complex(0,theta))
            cr.append(z)
        return rr, cr

    def computeRootsFromParams(self):
        self.realDenominatorRoots, self.complexDenominatorRoots = self.computeRootsFromParamsHelper(
            self.realDenominatorX, self.complexDenominatorX, self.complexDenominatorTheta, DiskRational.maxDenominatorRadius
        )

    @staticmethod
    def create(numRealNumeratorRoots, numComplexNumeratorRoots, numRealDenominatorRoots, numComplexDenominatorRoots):
        numerator = Util.filterPolynomialFromRoots(numRealNumeratorRoots, numComplexNumeratorRoots, DiskRational.maxNumeratorRadius)
        realDenominatorRoots, complexDenominatorRoots = Util.filterRoots(numRealDenominatorRoots, numComplexDenominatorRoots, DiskRational.maxDenominatorRadius)
        return DiskRational(numerator.coefficients, realDenominatorRoots, complexDenominatorRoots)

    def copy(self):
        r = DiskRational(self.numeratorCoef, self.realDenominatorRoots, self.complexDenominatorRoots)
        r.setFixedRoots(self.fixedNumeratorRoots, self.fixedDenominatorRoots)
        return r

    def incorporateRoots(self, numeratorRoots, denominatorRoots):
        p = Util.Polynomial(self.numeratorCoef)
        for r in numeratorRoots:
            p = p.mulMonomial(-r)
            p = p.mulMonomial(-r.conjugate())
        self.numeratorCoef = p.realifyCoef().coefficients
        rx, cx, ctheta = self.computeParamsFromRootsHelper([], denominatorRoots, self.maxDenominatorRadius)
        self.complexDenominatorX.extend(cx)
        self.complexDenominatorTheta.extend(ctheta)
        self.computeRootsFromParams()


    def evaluate(self, x):
        n = 0
        d = 1
        for i, a in enumerate(self.numeratorCoef):
            n += a * x**i
        for r in self.fixedNumeratorRoots:
            n *= (x - r)
            n *= (x - r.conjugate())
        for r in self.realDenominatorRoots:
            d *= (x - r)
        for r in self.complexDenominatorRoots:
            d *= (x - r)
            d *= (x - r.conjugate())
        for r in self.fixedDenominatorRoots:
            d *= (x - r)
            d *= (x - r.conjugate())
        return n / d

    '''def evaluatePartial(self, x):
        n = 1
        d = 1
        nParts = []
        dParts = []
        for r in self.realNumeratorRoots:
            t = (x - r)
            nParts.append(t)
            n *= t
        for r in self.complexNumeratorRoots:
            t = (x - r) * (x - r.conjugate())
            nParts.append(t)
            d *= t
        for r in self.realDenominatorRoots:
            t = (x - r)
            dParts.append(t)
            d *= t
        for r in self.complexDenominatorRoots:
            t = (x - r) * (x - r.conjugate())
            dParts.append(t)
            d *= t
        nParts = [n / t for t in nParts]
        dParts = [d / t for t in dParts]
        return nParts, dParts, n, d'''

    def evaluateDenominator(self, x):
        d = 1
        for r in self.realDenominatorRoots:
            d *= (x - r)
        for r in self.complexDenominatorRoots:
            d *= (x - r)
            d *= (x - r.conjugate())
        for r in self.fixedDenominatorRoots:
            d *= (x - r)
            d *= (x - r.conjugate())
        return d

    def updateCoef(self, deriv, rate):
        k = 0
        for i in range(len(self.numeratorCoef)):
            self.numeratorCoef[i] -= rate * deriv[k]
            k += 1
        for i in range(len(self.realDenominatorX)):
            self.realDenominatorX[i] -= rate * deriv[k]
            k += 1
        for i in range(len(self.complexDenominatorX)):
            self.complexDenominatorX[i] -= rate * deriv[k]
            self.complexDenominatorTheta[i] -= rate * deriv[k+1]
            t = self.complexDenominatorTheta[i] / (2*math.pi)
            self.complexDenominatorTheta[i] = (t - int(t)) * 2 * math.pi
            k += 2
        self.computeRootsFromParams()
        return self

    def getCoordinates(self):
        ret = []
        ret.extend(self.numeratorCoef)
        ret.extend(self.realDenominatorX)
        ret.extend(self.complexDenominatorX)
        ret.extend(self.complexDenominatorTheta)
        return ret

    def setFixedRoots(self, fixedNumeratorRoots, fixedDenominatorRoots):
        self.fixedNumeratorRoots = fixedNumeratorRoots
        self.fixedDenominatorRoots = fixedDenominatorRoots

    def evaluateFixedRoots(self, z):
        n = 1
        for r in self.fixedNumeratorRoots:
            n *= (z - r)
            n *= (z - r.conjugate())
        d = 1
        for r in self.fixedDenominatorRoots:
            d *= (z - r)
            d *= (z - r.conjugate())
        return n/d

    def derivative(self, grid, filter, sampleRate):
        deriv = [0 for i in range(len(self.numeratorCoef) +
                                  len(self.realDenominatorX) +
                                  len(self.complexDenominatorX) +
                                  len(self.complexDenominatorTheta))]
        for freq, filterValue in zip(grid, filter):
            z = cmath.exp(complex(0, -2 * math.pi * freq / sampleRate))
            y = self.evaluate(z)
            d = self.evaluateDenominator(z)
            f = abs(y)
            fixedPart = self.evaluateFixedRoots(z)
            prefactor = 2 * (f - filterValue) / f
            k = 0
            for i in range(len(self.numeratorCoef)):
                deriv[k] += prefactor * (y.conjugate() * fixedPart * z**i / d).real
                k += 1
            for i in range(len(self.realDenominatorX)):
                x = self.realDenominatorX[i]
                r = DiskRational.maxDenominatorRadius * (1 - math.exp(-x) / (1 + math.exp(-x)))
                z2 = z - r
                deriv[k] += prefactor * DiskRational.maxDenominatorRadius * (y * y / z2 * fixedPart * (-2*math.exp(-2*x)/(1+math.exp(-x))**2)).real
                k += 1
            for x, theta in zip(self.complexDenominatorX, self.complexDenominatorTheta):
                r = DiskRational.maxDenominatorRadius / (1 + math.exp(-x)) * cmath.exp(complex(0, theta))
                z2 = z - r
                z2rconj = z - r.conjugate()
                dX = DiskRational.maxDenominatorRadius * math.exp(-x)/(1+math.exp(-x))**2 * cmath.exp(complex(0, theta))
                dTheta = DiskRational.maxDenominatorRadius * complex(0, 1) * cmath.exp(complex(0, theta)) / (1+math.exp(-x))
                deriv[k] += prefactor * (y.conjugate() * fixedPart * y * (dX / z2 + dX.conjugate() / z2rconj)).real
                deriv[k+1] += prefactor * (y.conjugate() * fixedPart * y * (dTheta / z2 + dTheta.conjugate() / z2rconj)).real
                k += 2
        numGridPoints = len(grid)
        deriv = [t / numGridPoints for t in deriv]
        norm = functools.reduce(lambda n, t: n+t*t, deriv, 0)
        return deriv, math.sqrt(norm)

    def mutate(self, chanceOfMutation, mutationSize):
        newCoef = []
        for x in self.numeratorCoef:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
            newCoef.append(x)
        self.numeratorCoef = newCoef
        newCoef = []
        for x in self.realDenominatorX:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
            newCoef.append(x)
        self.realDenominatorX = newCoef
        newCoef = []
        for x in self.complexDenominatorX:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
            newCoef.append(x)
        self.complexDenominatorX = newCoef
        newCoef = []
        for x in self.complexDenominatorTheta:
            if random.uniform(0, 1) < chanceOfMutation:
                x *= 1 + random.uniform(-mutationSize/2, mutationSize/2)
            newCoef.append(x)
        self.complexDenominatorTheta = newCoef
        self.computeRootsFromParams()


    def crossover(self, other):
        newNumeratorCoef = []
        for x, y in zip(self.numeratorCoef, other.numeratorCoef):
            newNumeratorCoef.append(x if random.uniform(0, 1) < 0.5 else y)
        newRealDenominatorRoots = []
        for x, y in zip(self.realDenominatorRoots, other.realDenominatorRoots):
            newRealDenominatorRoots.append(x if random.uniform(0, 1) < 0.5 else y)
        newComplexDenominatorRoots = []
        for x, y in zip(self.complexDenominatorRoots, other.complexDenominatorRoots):
            newComplexDenominatorRoots.append(x if random.uniform(0, 1) < 0.5 else y)
        return DiskRational(newNumeratorCoef, newRealDenominatorRoots, newComplexDenominatorRoots)

if __name__ == '__main__':
    gridPoints = 200
    sampleRate = 96000

    grid = numpy.linspace(0, sampleRate/2, gridPoints)
    bspline = BSpline.BSpline(grid, BSpline.BSpline.getDefaultKnots())
    filter = bspline.getBasis(8)

    r = DiskRational.create(0, 2, 0, 2)
    r.setFixedRoots([complex(0.2, 0.3)], [])
    deriv, derivNorm = r.derivative(grid, filter, sampleRate)
    h = 0.001
    r2 = r.copy()
    r2.numeratorCoef[0] += h
    r2.computeRootsFromParams()
    d = (Util.objective(grid, filter, sampleRate, r2) - Util.objective(grid, filter, sampleRate, r)) / h
    print(len(r.numeratorCoef))
    print(len(r.realDenominatorX))
    print(len(r.complexDenominatorX))
    print(d)
    print(deriv[0])
    print(deriv)
