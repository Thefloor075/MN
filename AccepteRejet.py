
from math import *
import numpy as np
import matplotlib.pyplot as plt


def simps(f,a,b,N=50):
    '''Approximate the integral of f(x) from a to b by Simpson's rule.

    Simpson's rule approximates the integral \int_a^b f(x) dx by the sum:
    (dx/3) \sum_{k=1}^{N/2} (f(x_{2i-2} + 4f(x_{2i-1}) + f(x_{2i}))
    where x_i = a + i*dx and dx = (b - a)/N.

    Parameters
    ----------
    f : function
        Vectorized function of a single variable
    a , b : numbers
        Interval of integration [a,b]
    N : (even) integer
        Number of subintervals of [a,b]

    Returns
    -------
    float
        Approximation of the integral of f(x) from a to b using
        Simpson's rule with N subintervals of equal length.

    Examples
    --------
    >>> simps(lambda x : 3*x**2,0,1,10)
    1.0
    '''
    if N % 2 == 1:
        raise ValueError("N must be an even integer.")
    dx = (b-a)/N
    x = np.linspace(a,b,N+1)
    y = f(x)
    S = dx/3 * np.sum(y[0:-1:2] + 4*y[1::2] + y[2::2])
    return S

def integrale(f, a, b, n):
    somme = 0
    h = (b-a)/n
    x = a
    for _ in range(0, n + 1):
        somme += f(x)
        x += h
    return somme*h


def f(x):
	return (2 + np.sin(x)) * np.exp(- ( 2 + np.cos(3*x) + np.sin(2*x))*x)



X = np.linspace(-1, 10.0, num=10000)
F = f(X)

A = simps(f,0,200,10000000)
print(A)

plt.plot(X,F)
plt.show()
