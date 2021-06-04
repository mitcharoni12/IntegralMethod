import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from drawnow import *

def linFunction(x, a, b):
    return x**a + b

def expo(x):
    return x**1.8 + 90

for x in range(0, 100):
    xData[x] = x
    yData[x] = expo(x)

fitData, conv = curve_fit(linFunction, xData, yData)
a, b =  fitData
print(a)
print(b)