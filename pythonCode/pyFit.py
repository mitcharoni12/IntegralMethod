import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy import stats
from drawnow import *
import numpy as np
import math
from decimal import *

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

#def totalIntegral(x, Cs0, lambdaCs, Ba0, lambdaBa, La0, lambdaLa):
#    #Cs part of decay
#    f = Cs0 * (1 - np.exp(-lambdaCs * x))
#
#    #Ba part of decay
#    f += Ba0 * (1 - np.exp(-lambdaBa * x))
#
#    f += (((Cs0*lambdaBa)/(lambdaBa-lambdaCs))*(1 - np.exp(-lambdaCs * x)))
#    f += (((Cs0*lambdaCs)/(lambdaCs-lambdaBa))*(1 - np.exp(-lambdaBa * x)))
#
#    #La part of decay
#    f += La0 * ( 1.0 - np.exp(1 - lambdaLa * x ))
#
#    f += Ba0 * lambdaLa / (lambdaLa - lambdaBa ) * ( 1.0 - np.exp( - lambdaBa * x ))
#    f += Ba0 * lambdaBa / (lambdaBa - lambdaLa ) * ( 1.0 - np.exp( - lambdaLa * x ))
#
#    f += Cs0 * lambdaCs * lambdaBa / (lambdaCs - lambdaLa ) / ( lambdaBa - lambdaLa ) * ( 1.0 - np.exp( - lambdaLa * x ))
#    f += Cs0 * lambdaCs * lambdaLa / (lambdaCs - lambdaBa ) / ( lambdaLa - lambdaBa ) * ( 1.0 - np.exp( - lambdaBa * x ))
#    f += Cs0 * lambdaBa * lambdaLa / (lambdaBa - lambdaCs ) / ( lambdaLa - lambdaCs ) * ( 1.0 - np.exp( - lambdaCs * x ))
#
#    return f

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



integralDataFile = open("/home/rlmitchell43/integralMethod/IntegralMethod/multi/integralData.txt", "r")
eventYAxis = []
for x in range(100):
    binData = integralDataFile.readline()
    binData = float(binData)
    eventYAxis.append(binData)

integralBinDataFile = open("/home/rlmitchell43/integralMethod/IntegralMethod/multi/integralBinData.txt", "r")
binPoints = []
for x in range(100):
    binData = integralBinDataFile.readline()
    binData = float(binData)
    binPoints.append(binData)

def totalIntegral(params):
    Cs0 = params[0]
    lambdaCs = params[1]
    Ba0 = params[2]
    lambdaBa = params[3]
    La0 = params[4]
    lambdaLa = params[5]
    sd = params[6]

    #Cs part of decay
    f = Cs0 * (1.0 - np.exp(-lambdaCs * binPoints))

    #Ba part of decay
    f += Ba0 * (1.0 - np.exp(-lambdaBa * binPoints))

    f += (((Cs0*lambdaBa)/(lambdaBa-lambdaCs))*(1.0 - np.exp(-lambdaCs * binPoints)))
    f += (((Cs0*lambdaCs)/(lambdaCs-lambdaBa))*(1.0 - np.exp(-lambdaBa * binPoints)))

    #La part of decay
    f += La0 * ( 1.0 - np.exp(1.0 - lambdaLa * binPoints ))

    f += Ba0 * lambdaLa / (lambdaLa - lambdaBa ) * ( 1.0 - np.exp( - lambdaBa * binPoints ))
    f += Ba0 * lambdaBa / (lambdaBa - lambdaLa ) * ( 1.0 - np.exp( - lambdaLa * binPoints ))

    f += Cs0 * lambdaCs * lambdaBa / (lambdaCs - lambdaLa ) / ( lambdaBa - lambdaLa ) * ( 1.0 - np.exp( - lambdaLa * binPoints ))
    f += Cs0 * lambdaCs * lambdaLa / (lambdaCs - lambdaBa ) / ( lambdaLa - lambdaBa ) * ( 1.0 - np.exp( - lambdaBa * binPoints ))
    f += Cs0 * lambdaBa * lambdaLa / (lambdaBa - lambdaCs ) / ( lambdaLa - lambdaCs ) * ( 1.0 - np.exp( - lambdaCs * binPoints ))

    LL = -np.sum( stats.norm.logpdf(eventYAxis, loc=f, scale=sd ) )

    return (LL)

initParams = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

results = minimize(totalIntegral, initParams, method = 'Nelder-Mead')
print(results.x)

plt.plot(binPoints, eventYAxis)
plt.show()

boundArr = (
            (0, .01 * lambdaCs, 0, .01 * lambdaBa, 0, .01 * lambdaLa),
            (2 * Cs0, 100 * lambdaCs, 2 * Cs0, 100 * lambdaBa, 2 * Cs0, 100 * lambdaLa)
            )
'''
boundArr = (
            (0, .01 * lambdaCs, 0, .01 * lambdaBa),
            (2 * Cs0, 100 * lambdaCs, 2 * Cs0, 100 * lambdaBa)
            )
'''
#fitData, conv = curve_fit(LaDecayRegular, binCenters, n, bounds = boundArr, method = 'trf', maxfev = 10000000000000000000)
#Cs0Fit, lambdaCsFit, Ba0Fit, lambdaBaFit, La0Fit, lambdaLaFit = fitData
#Cs0Fit, lambdaCsFit, Ba0Fit, lambdaBaFit = fitData

#print("Cs0:           " + str(Cs0Fit))
#print("Cs Half Life:  " + str(np.log(2)/lambdaCsFit))
#print("Ba0:           " + str(Ba0Fit))
#print("Ba Half Life:  " + str(np.log(2)/lambdaBaFit))
#print("La0:           " + str(La0Fit))
#print("La Half Life:  " + str(np.log(2)/lambdaLaFit))

#plt.plot(binCenters, n)
#plt.show()

'''
x = np.linspace(0, 400, 10000)

plt.ylabel("Number Events")
plt.xlabel("Bin number")
plt.plot(x, totalRegular(x, Cs0, lambdaCs, Ba0, lambdaBa, La0, lambdaLa), label = "fit")
plt.plot(binXAxis, eventYAxis)
plt.show()
'''
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

