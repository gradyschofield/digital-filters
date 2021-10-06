import scipy.interpolate
import numpy
import matplotlib.pyplot as plt

class BSpline:

    @staticmethod
    def getDefaultKnots():
        return [0, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 48000]

    def getNumBases(self):
        return self.numBases

    def getBasis(self, i):
        return self.bases[i]

    def __init__(self, grid, knots, degree = 3):
        self.degree = degree
        self.knots = list(knots)
        for i in range(degree):
            self.knots.insert(0, self.knots[0])
            self.knots.append(self.knots[-1])
        self.numBases = len(self.knots)-degree-1
        self.bases = []
        for i in range(self.numBases):
            weights = [1 if i == j else 0 for j in range(self.numBases)]
            b1 = scipy.interpolate.BSpline(self.knots, weights, degree, extrapolate=False)
            x = numpy.linspace(self.knots[i], self.knots[i + degree + 1], 100)
            self.bases.append([b1(x) if x <= self.knots[i+degree+1] and x >= self.knots[i] else 0 for x in grid])


if __name__ == '__main__':
    degree = 3
    maxFreq = 48000
    grid = numpy.linspace(0, 48000, 1000)
    bspline = BSpline(grid, BSpline.getDefaultKnots(), 3)
    for i in range(bspline.getNumBases()):
        plt.plot(grid, bspline.getBasis(i))
    plt.show()
