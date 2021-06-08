import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from drawnow import *
import numpy as np
import math
from decimal import Decimal

Cs0, Ba0, La0 = 10000, 0, 0
lambdaCs, lambdaBa, lambdaLa = .69733, .060273, .016988

def CSDecayIntegral(x, Cs0, lambdaCs):
    f = Cs0 * (1 - np.exp(-lambdaCs * x))
    return f

def BADecayIntegral(x, Cs0, lambdaCs, Ba0, lambdaBa):
    f = Ba0 * (1 - np.exp(-lambdaBa * x))

    f += (((Cs0*lambdaBa)/(lambdaBa-lambdaCs))*(1 - np.exp(-lambdaCs * x)))
    f += (((Cs0*lambdaCs)/(lambdaCs-lambdaBa))*(1 - np.exp(-lambdaBa * x)))

    return f

def LADecayInegral(x, Cs0, lambdaCs, Ba0, lambdaBa, La0, lambdaLa):
    f = La0 * ( 1.0 - np.exp(1 - lambdaLa * x ))

    f += Ba0 * lambdaLa / (lambdaLa - lambdaBa ) * ( 1.0 - np.exp( - lambdaBa * x ))
    f += Ba0 * lambdaBa / (lambdaBa - lambdaLa ) * ( 1.0 - np.exp( - lambdaLa * x ))

    f += Cs0 * lambdaCs * lambdaBa / (lambdaCs - lambdaLa ) / ( lambdaBa - lambdaLa ) * ( 1.0 - np.exp( - lambdaLa * x ))
    f += Cs0 * lambdaCs * lambdaLa / (lambdaCs - lambdaBa ) / ( lambdaLa - lambdaBa ) * ( 1.0 - np.exp( - lambdaBa * x ))
    f += Cs0 * lambdaBa * lambdaLa / (lambdaBa - lambdaCs ) / ( lambdaLa - lambdaCs ) * ( 1.0 - np.exp( - lambdaCs * x ))

    return f

def totalIntegral(x, Cs0, lambdaCs, Ba0, lambdaBa, La0, lambdaLa):
    #Cs part of decay
    f = Cs0 * (1 - np.exp(-lambdaCs * x))

    #Ba part of decay
    f += Ba0 * (1 - np.exp(-lambdaBa * x))

    f += (((Cs0*lambdaBa)/(lambdaBa-lambdaCs))*(1 - np.exp(-lambdaCs * x)))
    f += (((Cs0*lambdaCs)/(lambdaCs-lambdaBa))*(1 - np.exp(-lambdaBa * x)))

    #La part of decay
    f += La0 * ( 1.0 - np.exp(1 - lambdaLa * x ))

    f += Ba0 * lambdaLa / (lambdaLa - lambdaBa ) * ( 1.0 - np.exp( - lambdaBa * x ))
    f += Ba0 * lambdaBa / (lambdaBa - lambdaLa ) * ( 1.0 - np.exp( - lambdaLa * x ))

    f += Cs0 * lambdaCs * lambdaBa / (lambdaCs - lambdaLa ) / ( lambdaBa - lambdaLa ) * ( 1.0 - np.exp( - lambdaLa * x ))
    f += Cs0 * lambdaCs * lambdaLa / (lambdaCs - lambdaBa ) / ( lambdaLa - lambdaBa ) * ( 1.0 - np.exp( - lambdaBa * x ))
    f += Cs0 * lambdaBa * lambdaLa / (lambdaBa - lambdaCs ) / ( lambdaLa - lambdaCs ) * ( 1.0 - np.exp( - lambdaCs * x ))

    return f

def CsDecayRegular(x, Cs0, lambdaCs):
    f = Cs0*lambdaCs*np.exp(- lambdaCs* x)
    return f

def BaDecayRegular(x, Cs0, lambdaCs, Ba0, lambdaBa):
    f = Ba0 * lambdaBa * np.exp(- lambdaBa*x)

    f += Cs0 * lambdaCs * lambdaBa * ((np.exp(- lambdaCs*x))/(lambdaBa-lambdaCs))
    f += Cs0 * lambdaCs * lambdaBa * ((np.exp(- lambdaBa*x))/(lambdaCs-lambdaBa))

    return f

def LaDecayRegular(x, Cs0, lambdaCs, Ba0, lambdaBa, La0, lambdaLa):
    f = (La0 * lambdaLa * (np.exp(-lambdaLa * x)))

    f += (Ba0 * lambdaBa *lambdaLa*((np.exp(- lambdaBa * x))/(lambdaLa - lambdaBa)))
    f += (Ba0 * lambdaBa *lambdaLa*((np.exp(- lambdaLa * x))/(lambdaBa - lambdaLa)))

    f += (Cs0 * lambdaCs * lambdaBa * lambdaLa*((np.exp(- lambdaCs * x))/((lambdaBa-lambdaCs)*(lambdaLa-lambdaCs))))
    f += (Cs0 * lambdaCs * lambdaBa * lambdaLa*((np.exp(- lambdaBa * x))/((lambdaCs-lambdaBa)*(lambdaLa-lambdaBa))))
    f += (Cs0 * lambdaCs * lambdaBa * lambdaLa*((np.exp(- lambdaLa * x))/((lambdaCs-lambdaLa)*(lambdaBa-lambdaLa))))

    return f

def totalRegular(x, Cs0, lambdaCs, Ba0, lambdaBa, La0, lambdaLa):
    #Cs part of decay
    f = Cs0*lambdaCs*np.exp(- lambdaCs* x)

    #Ba part of decay
    f += Ba0 * lambdaBa * np.exp(- lambdaBa*x)

    f += Cs0 * lambdaCs * lambdaBa * ((np.exp(- lambdaCs*x))/(lambdaBa-lambdaCs))
    f += Cs0 * lambdaCs * lambdaBa * ((np.exp(- lambdaBa*x))/(lambdaCs-lambdaBa))

    #La part of decay
    f += (La0 * lambdaLa * (np.exp(-lambdaLa * x)))

    f += (Ba0 * lambdaBa *lambdaLa*((np.exp(- lambdaBa * x))/(lambdaLa - lambdaBa)))
    f += (Ba0 * lambdaBa *lambdaLa*((np.exp(- lambdaLa * x))/(lambdaBa - lambdaLa)))

    f += (Cs0 * lambdaCs * lambdaBa * lambdaLa*((np.exp(- lambdaCs * x))/((lambdaBa-lambdaCs)*(lambdaLa-lambdaCs))))
    f += (Cs0 * lambdaCs * lambdaBa * lambdaLa*((np.exp(- lambdaBa * x))/((lambdaCs-lambdaBa)*(lambdaLa-lambdaBa))))
    f += (Cs0 * lambdaCs * lambdaBa * lambdaLa*((np.exp(- lambdaLa * x))/((lambdaCs-lambdaLa)*(lambdaBa-lambdaLa))))

    return f

myFile = open("/home/rlmitchell43/integralMethod/IntegralMethod/multi/regularValues.txt", "r")
binXAxis = []
eventYAxis = []
for x in range(1000):
    binXAxis.append(x)
    binData = myFile.readline()
    binData = Decimal(binData)
    eventYAxis.append(binData)

boundArr = (
            (0, .01 * lambdaCs, 0, .01 * lambdaBa, 0, .01 * lambdaLa),
            (2 * Cs0, 100 * lambdaCs, 2 * Cs0, 100 * lambdaBa, 2 * Cs0, 100 * lambdaLa)
            )
fitData, conv = curve_fit(totalRegular, binXAxis, binData, bounds = boundArr, method = 'trf')
Cs0Fit, lambdaCsFit, Ba0Fit, lambdaBaFit, La0Fit, lambdaLaFit = fitData

print(Cs0Fit)
print(lambdaCsFit)
print(Ba0Fit)
print(lambdaBaFit)
print(La0Fit)
print(lambdaLaFit)

x = np.linspace(0, 400, 10000)

plt.ylabel("Number Events")
plt.xlabel("Bin number")
plt.plot(x, totalRegular(x, Cs0, lambdaCs, Ba0, lambdaBa, La0, lambdaLa), label = "fit")
plt.plot(binXAxis, eventYAxis)
plt.show()

'''
xData = []
yData = []

for x in range(0, 100):
    y = expo(x)
    xData.append(x)
    yData.append(y)

fitData, conv = curve_fit(linFunction, xData, yData)
a, b =  fitData
print(a)
print(b)
'''

