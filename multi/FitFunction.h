#ifndef FITFUNCTION_H
#define FITFUNCTION_H

#include "TMath.h"

using namespace std;

class FitFunction{
public:
	typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);
	static Double_t CUDecayByActivity(Double_t *x, Double_t* par);
	static Double_t CUDecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t ZNDecayByActivity(Double_t *x, Double_t* par);
	static Double_t ZNDecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t GADecayByActivity(Double_t *x, Double_t* par);
	static Double_t GADecayByActivityIntegral(Double_t *x, Double_t* par);
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
	batemanFitFunctions[0] = CUDecayByActivity;
	integralFitFunctions[0] = CUDecayByActivityIntegral;
	batemanFitFunctions[1] = ZNDecayByActivity;
	integralFitFunctions[1] = ZNDecayByActivityIntegral;
	batemanFitFunctions[2] = GADecayByActivity;
	integralFitFunctions[2] = GADecayByActivityIntegral;
}

Double_t FitFunction::CUDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CU0 = par[0];
	Double_t lambdaCU = par[1];

	Double_t f = (CU0 * lambdaCU * (TMath::Exp(-lambdaCU * timeVar)));

	return f;
}

Double_t FitFunction::CUDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CU0 = par[0];
	Double_t lambdaCU = par[1];

	Double_t f = CU0 * (1.0 - TMath::Exp(-lambdaCU * timeVar));

	return f;
}

Double_t FitFunction::ZNDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CU0 = par[0];
	Double_t lambdaCU = par[1];
	Double_t ZN0 = par[2];
	Double_t lambdaZN = par[3];

	Double_t f = (ZN0 * lambdaZN * (TMath::Exp(-lambdaZN * timeVar)));

	f += (CU0 * lambdaZN * lambdaCU * ((TMath::Exp(-lambdaCU * timeVar)) / ((lambdaZN-lambdaCU))));
	f += (CU0 * lambdaZN * lambdaCU * ((TMath::Exp(-lambdaZN * timeVar)) / ((lambdaCU-lambdaZN))));

	return f;
}

Double_t FitFunction::ZNDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CU0 = par[0];
	Double_t lambdaCU = par[1];
	Double_t ZN0 = par[2];
	Double_t lambdaZN = par[3];

	Double_t f = ZN0 * (1.0 - TMath::Exp(-lambdaZN * timeVar));

	f += CU0 * lambdaZN * (1.0 - TMath::Exp(-lambdaCU * timeVar)) / ((lambdaZN-lambdaCU));
	f += CU0 * lambdaCU * (1.0 - TMath::Exp(-lambdaZN * timeVar)) / ((lambdaCU-lambdaZN));

	return f;
}

Double_t FitFunction::GADecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CU0 = par[0];
	Double_t lambdaCU = par[1];
	Double_t ZN0 = par[2];
	Double_t lambdaZN = par[3];
	Double_t GA0 = par[4];
	Double_t lambdaGA = par[5];

	Double_t f = (GA0 * lambdaGA * (TMath::Exp(-lambdaGA * timeVar)));

	f += (CU0 * lambdaGA * lambdaCU * lambdaZN * ((TMath::Exp(-lambdaCU * timeVar)) / ((lambdaZN-lambdaCU)*(lambdaGA-lambdaCU))));
	f += (CU0 * lambdaGA * lambdaCU * lambdaZN * ((TMath::Exp(-lambdaZN * timeVar)) / ((lambdaCU-lambdaZN)*(lambdaGA-lambdaZN))));
	f += (CU0 * lambdaGA * lambdaCU * lambdaZN * ((TMath::Exp(-lambdaGA * timeVar)) / ((lambdaCU-lambdaGA)*(lambdaZN-lambdaGA))));

	f += (ZN0 * lambdaGA * lambdaZN * ((TMath::Exp(-lambdaZN * timeVar)) / ((lambdaGA-lambdaZN))));
	f += (ZN0 * lambdaGA * lambdaZN * ((TMath::Exp(-lambdaGA * timeVar)) / ((lambdaZN-lambdaGA))));

	return f;
}

Double_t FitFunction::GADecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CU0 = par[0];
	Double_t lambdaCU = par[1];
	Double_t ZN0 = par[2];
	Double_t lambdaZN = par[3];
	Double_t GA0 = par[4];
	Double_t lambdaGA = par[5];

	Double_t f = GA0 * (1.0 - TMath::Exp(-lambdaGA * timeVar));

	f += CU0 * lambdaZN * lambdaGA * (1.0 - TMath::Exp(-lambdaCU * timeVar)) / ((lambdaZN-lambdaCU)*(lambdaGA-lambdaCU));
	f += CU0 * lambdaCU * lambdaGA * (1.0 - TMath::Exp(-lambdaZN * timeVar)) / ((lambdaCU-lambdaZN)*(lambdaGA-lambdaZN));
	f += CU0 * lambdaCU * lambdaZN * (1.0 - TMath::Exp(-lambdaGA * timeVar)) / ((lambdaCU-lambdaGA)*(lambdaZN-lambdaGA));

	f += ZN0 * lambdaGA * (1.0 - TMath::Exp(-lambdaZN * timeVar)) / ((lambdaGA-lambdaZN));
	f += ZN0 * lambdaZN * (1.0 - TMath::Exp(-lambdaGA * timeVar)) / ((lambdaZN-lambdaGA));

	return f;
}

#endif