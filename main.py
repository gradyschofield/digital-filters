# This is a sample Python script.
import random
import math
import cmath
import matplotlib.pyplot as plt
import Util
# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.

class AllPassFilter:
    denominatorParams = []
    sampleRate = 0

    def __init__(self, coef, sampleRate):
        self.denominatorParams = coef
        self.sampleRate = 96000.0

    def eval(self, f):
        t = 2 * math.pi * f / self.sampleRate
        z = complex(math.cos(t), math.sin(t))
        d = 0
        n = 0
        for i in range(len(self.denominatorParams)):
            idx = len(self.denominatorParams) - 1 - i
            d += self.denominatorParams[idx] * z**-i
            n += self.denominatorParams[i] * z**-i
        return n / d

    # Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    p = Util.filterPolynomialFromRoots(2, 10, 0.99999)
    filter = AllPassFilter(p.coefficients, 96000.0)
    x = []
    y = []
    gd = []
    lastP = None
    for i in range(100, 48000, 100):
        x.append(i)
        p = cmath.phase(filter.eval(i))
        if lastP:
            while p - lastP > math.pi:
                p -= 2 * math.pi
            while p - lastP < -math.pi:
                p += 2 * math.pi
            #gd.append((p - lastP) / (2 * math.pi / 96000))
            gd.append(180 * (p - lastP) / math.pi)
            #sin( 2 pi f t + phaseShift)
        y.append(180 / math.pi * p)
        lastP = p
    #plt.plot(x, y)
    plt.plot(x[1:len(x)], gd)
    plt.show()


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
