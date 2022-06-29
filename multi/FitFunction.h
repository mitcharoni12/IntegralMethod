#ifndef FITFUNCTION_H
#define FITFUNCTION_H

#include "TMath.h"

using namespace std;

class FitFunction{
public:
	typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);
	static Double_t Test1DecayByActivity(Double_t *x, Double_t* par);
	static Double_t Test1DecayByActivityIntegral(Double_t *x, Double_t* par);
	FitFunction(Int_t numElements);
	decayFunction* GetBatemanFitFunctions(){return batemanFitFunctions;}
	decayFunction* GetIntegralFitFunctions(){return integralFitFunctions;}
private:
	Int_t numElements;
	decayFunction* batemanFitFunctions;
	decayFunction* integralFitFunctions;
};

FitFunction::FitFunction(Int_t numElements)
{
	this->numElements = numElements;
	batemanFitFunctions = new decayFunction [numElements];
	integralFitFunctions = new decayFunction [numElements];
	batemanFitFunctions[0] = Test1DecayByActivity;
	integralFitFunctions[0] = Test1DecayByActivityIntegral;
}

Double_t FitFunction::Test1DecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test10 = par[0];
	Double_t lambdaTest1 = par[1];

	Double_t f = (Test10 * lambdaTest1 * (TMath::Exp(-lambdaTest1 * timeVar)));

	return f;
}

Double_t FitFunction::Test1DecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test10 = par[0];
	Double_t lambdaTest1 = par[1];

	Double_t f = Test10 * (1.0 - TMath::Exp(-lambdaTest1 * timeVar));

	return f;
}

#endif