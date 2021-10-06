import random
import math
import cmath
import numpy

def randomComplex(maxRadius):
    r = random.uniform(0, maxRadius)
    angle = random.uniform(0, 2 * math.pi - 0.00001)
    return complex(r * math.cos(angle), r * math.sin(angle))

def randomComplexHalfAnnulus(innerRadius, outerRadius):
    r = random.uniform(innerRadius, outerRadius)
    angle = random.uniform(1E-5, math.pi - 1E-5)
    return complex(r * math.cos(angle), r * math.sin(angle))

class Polynomial:
    coefficients = []

    def mulMonomial(self, a):
        newCoef = list(self.coefficients)
        newCoef.insert(0, 0)
        for i in range(0, len(newCoef)-1):
            newCoef[i] += a * newCoef[i+1]
        return Polynomial(newCoef)

    @staticmethod
    def fromRoots(roots):
        ret = Polynomial([1])
        for r in roots:
            ret = ret.mulMonomial(-r)
        return ret

    def realifyCoef(self):
        newCoef = []
        for x in self.coefficients:
            if isinstance(x, complex):
                newCoef.append(x.real)
            else:
                newCoef.append(float(x))
        return Polynomial(newCoef)

    def __init__(self, coef):
        self.coefficients = coef

def filterPolynomialFromRoots(numReal, numComplex, radius):
    roots = []
    for i in range(numReal):
        roots.append(complex(random.uniform(-radius, radius), 0))
    for i in range(numComplex):
        r = randomComplex(radius)
        roots.append(r);
        roots.append(r.conjugate())
    #for r in roots:
    #    print(r)
    return Polynomial.fromRoots(roots).realifyCoef()

def filterRoots(numReal, numComplex, maxRadius):
    realRoots = []
    complexRoots = []
    for i in range(numReal):
        realRoots.append(random.uniform(-maxRadius, maxRadius))
    for i in range(numComplex):
        complexRoots.append(randomComplex(maxRadius))
    return realRoots, complexRoots

def filterRootsHalfAnnulus(numComplex, innerRadius, outerRadius):
    roots = []
    for i in range(numComplex):
        roots.append(randomComplexHalfAnnulus(innerRadius, outerRadius))
    return roots

def filterPolynomialWithCircularRootPattern(numRoots, radius, phase, boost):
    roots = []
    start = math.pi / 480000 * 10000
    stop = math.pi / 480000 * 20000
    grid, spacing = numpy.linspace(0, math.pi, numRoots, endpoint=False, retstep=True)
    grid = [x + spacing * phase for x in grid[1:len(grid)]]
    for x in grid:
        if x > start and x < stop:
            r = boost * cmath.exp(complex(0,x))
        else:
            r = radius * cmath.exp(complex(0, x))
        roots.append(r)
        roots.append(r.conjugate())
    return Polynomial.fromRoots(roots).realifyCoef()

def objective(grid, filter, sampleRate, rational):
    obj = 0
    sampleRateInv = 1.0 / (sampleRate+2)
    numGridPoints = len(grid)
    for freq, filterValue in zip(grid, filter):
        normalizedCoord = 2 * freq * sampleRateInv - 1 + sampleRateInv
        weight = 1#1 / math.sqrt(1 - normalizedCoord**2)
        x = cmath.exp(complex(0, -2 * math.pi * freq / sampleRate))
        y = rational.evaluate(x)
        f = abs(y)
        obj += weight * (f - filterValue) ** 2
    return obj / numGridPoints

