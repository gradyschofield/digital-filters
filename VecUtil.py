import math

def dot(v1, v2):
    if len(v1) != len(v2):
        raise "dot called on vectors not of the same length"
    d = 0
    for x1, x2 in zip(v1, v2):
        d += x1 * x2
    return d

def norm(v):
    n = 0
    for x in v:
        n += x**2
    return math.sqrt(n)

def maxLength(v, len, scale = 1.0):
    if len == 0:
        return v
    n = norm(v)
    if n * scale > len:
        return [len*x/n for x in v]
    return v

