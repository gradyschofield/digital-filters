import math
import numpy
import matplotlib.pyplot as plt

class Chebyshev:
    coefficients = []

    def __init__(self, coefficients):
        self.coefficients = coefficients

    @staticmethod
    def getZeros(degree):
        roots = []
        for i in range(0, degree):
            roots.append(math.cos(math.pi * (i + 0.5) / degree))
        return roots

    @staticmethod
    def getIntegrationWeight(degree):
        return math.pi / degree

    @staticmethod
    def projectStep(start, stop, min, max, gain, isFilter, degree):
        normalizedStart = 2 * (start - min) / (max-min) - 1
        normalizedStop = 2 * (stop - min) / (max-min) - 1
        step = lambda x: gain if x >= normalizedStart and x <= normalizedStop else 0 if isFilter else 1
        zeros = Chebyshev.getZeros(degree+1)
        weight = Chebyshev.getIntegrationWeight(degree+1)
        coefficients = []
        aK = math.pi / (degree + 2)
        for deg in range(degree+1):
            c = Chebyshev([0 if i < deg else 1 for i in range(deg + 1)])
            integral = 0
            for x in zeros:
                integral += weight * c.evaluate(x) * step(x)
            integral *= 1 / math.pi if deg == 0 else 2 / math.pi
            i = deg
            gI = ((1-i/float(degree+2)) * math.sin(aK)*math.cos(i*aK) + 1/(degree+2)*math.cos(aK)*math.sin(i*aK))/math.sin(aK)
            coefficients.append(gI * integral)
        return Chebyshev(coefficients)

    def evaluate(self, x):
        t1 = 1
        if len(self.coefficients) == 1:
            return self.coefficients[0] * t1
        t2 = x
        ret = self.coefficients[1] * t2 + self.coefficients[0] * t1
        for i in range(2,len(self.coefficients)):
            tmp = t2
            t2 = 2 * x * t2 - t1
            t1 = tmp
            ret += self.coefficients[i] * t2
        return ret

    def getRoots(self):
        coef = [0 for i in range(len(self.coefficients))]
        coef[0] = self.coefficients[0]
        coef[1] = self.coefficients[1]
        coef1 = [1, 0]
        if len(self.coefficients) > 1:
            coef2 = [0, 1]

            for i in range(2, len(self.coefficients)):
                coef1.append(0)
                coef2.insert(0, 0)
                t = [2*t2 - t1 for t2, t1, in zip(coef2, coef1)]
                coef1 = coef2
                coef2 = t
                for i in range(len(t)):
                    coef[i] += self.coefficients[i] * t[i]
        coef.reverse()
        return numpy.roots(coef)

    def evaluateOnScale(self, globalX, globalMin, globalMax):
        y = []
        for x in globalX:
            xNorm = 2 * (x - globalMin)/(globalMax - globalMin) - 1
            y.append(self.evaluate(xNorm))
        return y

    def plotData(self, min, max):
        x = []
        y = []
        steps = 200
        inc = (max-min)/steps
        for i in range(0, steps):
            t = inc * i + min
            x.append(t)
            y.append(self.evaluate(t))
        return x, y

    def plotDataOnScale(self, min, max, globalMin, globalMax):
        x = []
        y = []
        steps = 200
        max = 2 * (max - globalMin) / (globalMax - globalMin) - 1
        min = 2 * (min - globalMin) / (globalMax - globalMin) - 1
        inc = (max-min)/steps
        for i in range(0, steps):
            t = inc * i + min
            x.append(((t + 1) / 2) * (globalMax - globalMin) + globalMin)
            y.append(self.evaluate(t))
        return x, y

if __name__ == '__main__':
    c = Chebyshev.projectStep(500, 1000, 0, 20000, 6, False, 180)
    x, y = c.plotDataOnScale(0, 2000, 0, 20000)
    plt.plot(x, y)
    '''for deg in range(0,5):
        zeros = Chebyshev.getZeros(deg+1)
        weight = Chebyshev.getIntegrationWeight(deg+1)
        coef = [0 if i < deg else 1 for i in range(deg+1)]
        c = Chebyshev(coef)
        integral = 0
        for x in zeros:
            integral += c.evaluate(x)**2 * weight
        print(integral)
        x = []
        y = []
        inc = 0.01
        for i in range(0, int(2 / inc)):
            t = inc * i - 1
            x.append(t)
            y.append(c.evaluate(t))
        plt.plot(x, y)
        '''
    plt.show()

