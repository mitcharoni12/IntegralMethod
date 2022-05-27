#ifndef FITFUNCTION_H
#define FITFUNCTION_H

#include "TMath.h"

using namespace std;

class FitFunction{
public:
	typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);
	static Double_t TestDecayByActivity(Double_t *x, Double_t* par);
	static Double_t TestDecayByActivityIntegral(Double_t *x, Double_t* par);
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
	batemanFitFunctions[0] = TestDecayByActivity;
	integralFitFunctions[0] = TestDecayByActivityIntegral;
}

Double_t FitFunction::TestDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test0 = par[0];
	Double_t lambdaTest = par[1];

	Double_t f = (Test0 * lambdaTest * (TMath::Exp(-lambdaTest * timeVar)));

	return f;
}

Double_t FitFunction::TestDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test0 = par[0];
	Double_t lambdaTest = par[1];

	Double_t f = Test0 * (1.0 - TMath::Exp(-lambdaTest * timeVar));

	return f;
}

#endif