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
	static Double_t Test3DecayByActivity(Double_t *x, Double_t* par);
	static Double_t Test3DecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t Test4DecayByActivity(Double_t *x, Double_t* par);
	static Double_t Test4DecayByActivityIntegral(Double_t *x, Double_t* par);
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
	batemanFitFunctions[2] = Test3DecayByActivity;
	integralFitFunctions[2] = Test3DecayByActivityIntegral;
	batemanFitFunctions[3] = Test4DecayByActivity;
	integralFitFunctions[3] = Test4DecayByActivityIntegral;
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

Double_t FitFunction::Test3DecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test10 = par[0];
	Double_t lambdaTest1 = par[1];
	Double_t Test20 = par[2];
	Double_t lambdaTest2 = par[3];
	Double_t Test30 = par[4];
	Double_t lambdaTest3 = par[5];

	Double_t f = (Test30 * lambdaTest3 * (TMath::Exp(-lambdaTest3 * timeVar)));

	f += (Test10 * lambdaTest3 * lambdaTest1 * lambdaTest2 * ((TMath::Exp(-lambdaTest1 * timeVar)) / ((lambdaTest2-lambdaTest1)*(lambdaTest3-lambdaTest1))));
	f += (Test10 * lambdaTest3 * lambdaTest1 * lambdaTest2 * ((TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest1-lambdaTest2)*(lambdaTest3-lambdaTest2))));
	f += (Test10 * lambdaTest3 * lambdaTest1 * lambdaTest2 * ((TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest1-lambdaTest3)*(lambdaTest2-lambdaTest3))));

	f += (Test20 * lambdaTest3 * lambdaTest2 * ((TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest3-lambdaTest2))));
	f += (Test20 * lambdaTest3 * lambdaTest2 * ((TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest2-lambdaTest3))));

	return f;
}

Double_t FitFunction::Test3DecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test10 = par[0];
	Double_t lambdaTest1 = par[1];
	Double_t Test20 = par[2];
	Double_t lambdaTest2 = par[3];
	Double_t Test30 = par[4];
	Double_t lambdaTest3 = par[5];

	Double_t f = Test30 * (1.0 - TMath::Exp(-lambdaTest3 * timeVar));

	f += Test10 * lambdaTest2 * lambdaTest3 * (1.0 - TMath::Exp(-lambdaTest1 * timeVar)) / ((lambdaTest2-lambdaTest1)*(lambdaTest3-lambdaTest1));
	f += Test10 * lambdaTest1 * lambdaTest3 * (1.0 - TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest1-lambdaTest2)*(lambdaTest3-lambdaTest2));
	f += Test10 * lambdaTest1 * lambdaTest2 * (1.0 - TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest1-lambdaTest3)*(lambdaTest2-lambdaTest3));

	f += Test20 * lambdaTest3 * (1.0 - TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest3-lambdaTest2));
	f += Test20 * lambdaTest2 * (1.0 - TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest2-lambdaTest3));

	return f;
}

Double_t FitFunction::Test4DecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test10 = par[0];
	Double_t lambdaTest1 = par[1];
	Double_t Test20 = par[2];
	Double_t lambdaTest2 = par[3];
	Double_t Test30 = par[4];
	Double_t lambdaTest3 = par[5];
	Double_t Test40 = par[6];
	Double_t lambdaTest4 = par[7];

	Double_t f = (Test40 * lambdaTest4 * (TMath::Exp(-lambdaTest4 * timeVar)));

	f += (Test10 * lambdaTest4 * lambdaTest1 * lambdaTest2 * lambdaTest3 * ((TMath::Exp(-lambdaTest1 * timeVar)) / ((lambdaTest2-lambdaTest1)*(lambdaTest3-lambdaTest1)*(lambdaTest4-lambdaTest1))));
	f += (Test10 * lambdaTest4 * lambdaTest1 * lambdaTest2 * lambdaTest3 * ((TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest1-lambdaTest2)*(lambdaTest3-lambdaTest2)*(lambdaTest4-lambdaTest2))));
	f += (Test10 * lambdaTest4 * lambdaTest1 * lambdaTest2 * lambdaTest3 * ((TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest1-lambdaTest3)*(lambdaTest2-lambdaTest3)*(lambdaTest4-lambdaTest3))));
	f += (Test10 * lambdaTest4 * lambdaTest1 * lambdaTest2 * lambdaTest3 * ((TMath::Exp(-lambdaTest4 * timeVar)) / ((lambdaTest1-lambdaTest4)*(lambdaTest2-lambdaTest4)*(lambdaTest3-lambdaTest4))));

	f += (Test20 * lambdaTest4 * lambdaTest2 * lambdaTest3 * ((TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest3-lambdaTest2)*(lambdaTest4-lambdaTest2))));
	f += (Test20 * lambdaTest4 * lambdaTest2 * lambdaTest3 * ((TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest2-lambdaTest3)*(lambdaTest4-lambdaTest3))));
	f += (Test20 * lambdaTest4 * lambdaTest2 * lambdaTest3 * ((TMath::Exp(-lambdaTest4 * timeVar)) / ((lambdaTest2-lambdaTest4)*(lambdaTest3-lambdaTest4))));

	f += (Test30 * lambdaTest4 * lambdaTest3 * ((TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest4-lambdaTest3))));
	f += (Test30 * lambdaTest4 * lambdaTest3 * ((TMath::Exp(-lambdaTest4 * timeVar)) / ((lambdaTest3-lambdaTest4))));

	return f;
}

Double_t FitFunction::Test4DecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t Test10 = par[0];
	Double_t lambdaTest1 = par[1];
	Double_t Test20 = par[2];
	Double_t lambdaTest2 = par[3];
	Double_t Test30 = par[4];
	Double_t lambdaTest3 = par[5];
	Double_t Test40 = par[6];
	Double_t lambdaTest4 = par[7];

	Double_t f = Test40 * (1.0 - TMath::Exp(-lambdaTest4 * timeVar));

	f += Test10 * lambdaTest2 * lambdaTest3 * lambdaTest4 * (1.0 - TMath::Exp(-lambdaTest1 * timeVar)) / ((lambdaTest2-lambdaTest1)*(lambdaTest3-lambdaTest1)*(lambdaTest4-lambdaTest1));
	f += Test10 * lambdaTest1 * lambdaTest3 * lambdaTest4 * (1.0 - TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest1-lambdaTest2)*(lambdaTest3-lambdaTest2)*(lambdaTest4-lambdaTest2));
	f += Test10 * lambdaTest1 * lambdaTest2 * lambdaTest4 * (1.0 - TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest1-lambdaTest3)*(lambdaTest2-lambdaTest3)*(lambdaTest4-lambdaTest3));
	f += Test10 * lambdaTest1 * lambdaTest2 * lambdaTest3 * (1.0 - TMath::Exp(-lambdaTest4 * timeVar)) / ((lambdaTest1-lambdaTest4)*(lambdaTest2-lambdaTest4)*(lambdaTest3-lambdaTest4));

	f += Test20 * lambdaTest3 * lambdaTest4 * (1.0 - TMath::Exp(-lambdaTest2 * timeVar)) / ((lambdaTest3-lambdaTest2)*(lambdaTest4-lambdaTest2));
	f += Test20 * lambdaTest2 * lambdaTest4 * (1.0 - TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest2-lambdaTest3)*(lambdaTest4-lambdaTest3));
	f += Test20 * lambdaTest2 * lambdaTest3 * (1.0 - TMath::Exp(-lambdaTest4 * timeVar)) / ((lambdaTest2-lambdaTest4)*(lambdaTest3-lambdaTest4));

	f += Test30 * lambdaTest4 * (1.0 - TMath::Exp(-lambdaTest3 * timeVar)) / ((lambdaTest4-lambdaTest3));
	f += Test30 * lambdaTest3 * (1.0 - TMath::Exp(-lambdaTest4 * timeVar)) / ((lambdaTest3-lambdaTest4));

	return f;
}

#endif