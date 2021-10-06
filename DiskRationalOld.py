import cmath
import math
import functools

import Util
import BaseRational

class DiskRational(BaseRational.BaseRational):
    maxNumeratorRadius = 1.2
    maxDenominatorRadius = 0.999
    realNumeratorRoots = []
    complexNumeratorRoots = []
    realDenominatorRoots = []
    complexDenominatorRoots = []

    realNumeratorX = []
    complexNumeratorX = []
    complexNumeratorTheta = []
    realDenominatorX = []
    complexDenominatorX = []
    complexDenominatorTheta = []

    def __init__(self, realNumeratorRoots, complexNumeratorRoots, realDenominatorRoots, complexDenominatorRoots):
        self.realNumeratorRoots = list(realNumeratorRoots)
        self.complexNumeratorRoots = list(complexNumeratorRoots)
        self.realDenominatorRoots = list(realDenominatorRoots)
        self.complexDenominatorRoots = list(complexDenominatorRoots)
        self.computeParamsFromRoots()

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
        self.realNumeratorX, \
            self.complexNumeratorX, \
            self.complexNumeratorTheta = self.computeParamsFromRootsHelper(self.realNumeratorRoots, self.complexNumeratorRoots, DiskRational.maxNumeratorRadius)
        self.realDenominatorX, \
            self.complexDenominatorX, \
            self.complexDenominatorTheta = self.computeParamsFromRootsHelper(self.realDenominatorRoots, self.complexDenominatorRoots, DiskRational.maxDenominatorRadius)

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
        self.realNumeratorRoots, self.complexNumeratorRoots = self.computeRootsFromParamsHelper(
            self.realNumeratorX, self.complexNumeratorX, self.complexNumeratorTheta, DiskRational.maxNumeratorRadius
        )
        self.realDenominatorRoots, self.complexDenominatorRoots = self.computeRootsFromParamsHelper(
            self.realDenominatorX, self.complexDenominatorX, self.complexDenominatorTheta, DiskRational.maxDenominatorRadius
        )

    @staticmethod
    def create(numRealNumeratorRoots, numComplexNumeratorRoots, numRealDenominatorRoots, numComplexDenominatorRoots):
        realNumeratorRoots, complexNumeratorRoots = Util.filterRoots(numRealNumeratorRoots, numComplexNumeratorRoots, DiskRational.maxNumeratorRadius)
        realDenominatorRoots, complexDenominatorRoots = Util.filterRoots(numRealDenominatorRoots, numComplexDenominatorRoots, DiskRational.maxDenominatorRadius)
        return DiskRational(realNumeratorRoots, complexNumeratorRoots, realDenominatorRoots, complexDenominatorRoots)

    def copy(self):
        return DiskRational(self.realNumeratorRoots, self.complexNumeratorRoots, self.realDenominatorRoots, self.complexDenominatorRoots)

    def evaluate(self, x):
        n = 1
        d = 1
        for r in self.realNumeratorRoots:
            n *= (x - r)
        for r in self.complexNumeratorRoots:
            n *= (x - r)
            n *= (x - r.conjugate())
        for r in self.realDenominatorRoots:
            d *= (x - r)
        for r in self.complexDenominatorRoots:
            d *= (x - r)
            d *= (x - r.conjugate())
        return n / d

    def evaluatePartial(self, x):
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
        return nParts, dParts, n, d

    def evaluateDenominator(self, x):
        d = 1
        for r in self.realDenominatorRoots:
            d *= (x - r)
        for r in self.complexDenominatorRoots:
            d *= (x - r)
            d *= (x - r.conjugate())
        return d

    def updateCoef(self, deriv, rate):
        k = 0
        for i in range(len(self.realNumeratorX)):
            self.realNumeratorX[i] -= rate * deriv[k]
            k += 1
        for i in range(len(self.complexNumeratorX)):
            self.complexNumeratorX[i] -= rate * deriv[k]
            self.complexNumeratorTheta[i] -= rate * deriv[k+1]
            t = self.complexNumeratorTheta[i] / (2*math.pi)
            self.complexNumeratorTheta[i] = (t - int(t)) * 2 * math.pi
            k += 2
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

    def derivative(self, grid, filter, sampleRate):
        deriv = [0 for i in range(len(self.realNumeratorX) +
                                  len(self.complexNumeratorX) +
                                  len(self.complexNumeratorTheta) +
                                  len(self.realDenominatorX) +
                                  len(self.complexDenominatorX) +
                                  len(self.complexDenominatorTheta))]
        for freq, filterValue in zip(grid, filter):
            z = cmath.exp(complex(0, -2 * math.pi * freq / sampleRate))
            nParts, dParts, n, d = self.evaluatePartial(z)
            y = n / d
            f = abs(y)
            prefactor = 2 * (f - filterValue) / f
            k = 0
            param = 0
            for i in range(len(self.realNumeratorX)):
                x = self.realNumeratorX[i]
                deriv[k] += prefactor * (y * nParts[param] / d * (-2*math.exp(-2*x)/(1+math.exp(-x))**2)).real
                k += 1
                param += 1
            for x, theta in zip(self.complexNumeratorX, self.complexNumeratorTheta):
                r = DiskRational.maxNumeratorRadius / (1 + math.exp(-x)) * cmath.exp(complex(0, theta))
                z2 = z - r
                dX = -math.exp(-x)/(1+math.exp(-x))**2 * cmath.exp(complex(0, theta))
                dTheta = complex(0, 1) * cmath.exp(complex(0, theta)) / (1+math.exp(-x))
                deriv[k] += prefactor * (y * nParts[param] / d * 2 * (z2 * dX).real).real
                deriv[k+1] += prefactor * (y * nParts[param] / d * 2 * (z2 * dTheta).real).real
                k += 2
                param += 1
            for i in range(len(self.realDenominatorX)):
                x = self.realDenominatorX[i]
                r = DiskRational.maxDenominatorRadius * (1 - math.exp(-x) / (1 + math.exp(-x)))
                z2 = z - r
                deriv[k] += prefactor * (y * n / (d * z2) * (-2*math.exp(-2*x)/(1+math.exp(-x))**2)).real
                k += 1
            for x, theta in zip(self.complexDenominatorX, self.complexDenominatorTheta):
                r = DiskRational.maxDenominatorRadius / (1 + math.exp(-x)) * cmath.exp(complex(0, theta))
                z2 = z - r
                dX = -math.exp(-x)/(1+math.exp(-x))**2 * cmath.exp(complex(0, theta))
                dTheta = complex(0, 1) * cmath.exp(complex(0, theta)) / (1+math.exp(-x))
                deriv[k] += prefactor * (y * n / d * 2 * (dX / z2).real).real
                deriv[k+1] += prefactor * (y * n / d * 2 * (dTheta / z2).real).real
                k += 2
        numGridPoints = len(grid)
        deriv = [t / numGridPoints for t in deriv]
        norm = functools.reduce(lambda n, t: n+t*t, deriv, 0)
        return deriv, math.sqrt(norm)