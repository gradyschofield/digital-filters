import cmath
import math

class BaseRational:
    def plotData(self, min, max):
        x = []
        y = []
        steps = 200
        inc = (max - min) / steps
        for i in range(0, steps):
            t = inc * i + min
            x.append(t)
            y.append(self.evaluate(t))
        return x, y

    def plotData(self, min, max, sampleRate, steps = 200):
        x = []
        y = []
        #min = min / (2 * math.pi) * sampleRate
        #max = max / (2 * math.pi) * sampleRate
        inc = (max-min)/steps
        for i in range(0, steps):
            t = inc * i + min
            z = cmath.exp(complex(0,-2 * math.pi / sampleRate * t))
            x.append(t)
            y.append(abs(self.evaluate(z)))
        return x, y

