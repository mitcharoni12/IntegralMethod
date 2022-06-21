#ifndef FITFUNCTION_H
#define FITFUNCTION_H

#include "TMath.h"

using namespace std;

class FitFunction{
public:
	typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);
	static Double_t Test1DecayByActivity(Double_t *x, Double_t* par);
	static Double_t Test1DecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t Test2DecayByActivity(Double_t *x, Double_t* par);
	static Double_t Test2DecayByActivityIntegral(Double_t *x, Double_t* par);
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
	batemanFitFunctions[1] = Test2DecayByActivity;
	integralFitFunctions[1] = Test2DecayByActivityIntegral;
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

Double_t FitFunction::Test2DecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test10 = par[0];
	Double_t lambdaTest1 = par[1];
	Double_t Test20 = par[2];
	Double_t lambdaTest2 = par[3];

	Double_t f = (Test20 * lambdaTest2 * (TMath::Exp(-lambdaTest2 * timeVar)));

	f += (Test10 * lambdaTest2 * lambdaTest1 * ((TMath::Exp(-lambdaTest1 * timeVar)) / ((lambdaTest2-lambdaTest1))));
	f += (Test10 * lambdaTest2 * lambdaTest1 * ((TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest1-lambdaTest2))));

	return f;
}

Double_t FitFunction::Test2DecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test10 = par[0];
	Double_t lambdaTest1 = par[1];
	Double_t Test20 = par[2];
	Double_t lambdaTest2 = par[3];

	Double_t f = Test20 * (1.0 - TMath::Exp(-lambdaTest2 * timeVar));

	f += Test10 * lambdaTest2 * (1.0 - TMath::Exp(-lambdaTest1 * timeVar)) / ((lambdaTest2-lambdaTest1));
	f += Test10 * lambdaTest1 * (1.0 - TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest1-lambdaTest2));

	return f;
}

#endif