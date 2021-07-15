#ifndef FITFUNCTION_H
#define FITFUNCTION_H

#include "TMath.h"

typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);

using namespace std;

class FitFunction{
private:
	Int_t numElements;
	decayFunction* batemanFitFunctions;
	decayFunction* integralFitFunctions;
	decayFunction* fitFunctions;
public:
	Double_t CSDecayByActivity(Double_t *x, Double_t par);
	Double_t CSDecayByActivityIntegral(Double_t *x, Double_t par);
	Double_t BADecayByActivity(Double_t *x, Double_t par);
	Double_t BADecayByActivityIntegral(Double_t *x, Double_t par);
	Double_t LADecayByActivity(Double_t *x, Double_t par);
	Double_t LADecayByActivityIntegral(Double_t *x, Double_t par);
};

FitFunction::FitFunction(Int_t numElements)
{
	this->numElements = numElements;
	batemanFitFunctions = new decayFunction [numElements];
	integralFitFunctions = new decayFunction [numElements];
	fitFunctions = new decayFunction [numElements*2];
	fitFunction[0] = CSDecayByActivity;
	fitFunction[1] = CSDecayByActivityIntegral;
	batemanFitFunctions[0] = CSDecayByActivity;
	integralFitFunction[0] = CSDecayByActivityIntegral;
	fitFunction[2] = BADecayByActivity;
	fitFunction[3] = BADecayByActivityIntegral;
	batemanFitFunctions[1] = BADecayByActivity;
	integralFitFunction[1] = BADecayByActivityIntegral;
	fitFunction[4] = LADecayByActivity;
	fitFunction[5] = LADecayByActivityIntegral;
	batemanFitFunctions[2] = LADecayByActivity;
	integralFitFunction[2] = LADecayByActivityIntegral;
}

Double_t FitFunction::CSDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CS0 = par[0];
	Double_t lambdaCS = par[1];

	Double_t f = (CS0 * lambdaCS * (TMath::Exp(-lambdaCS * timeVar)));

	return f;
}

Double_t FitFunction::CSDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CS0 = par[0];
	Double_t lambdaCS = par[1];

	Double_t f = CS0 * (1.0 - TMath::Exp(-lambdaCS * timeVar));

	return f;
}

Double_t FitFunction::BADecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CS0 = par[0];
	Double_t lambdaCS = par[1];
	Double_t BA0 = par[2];
	Double_t lambdaBA = par[3];

	Double_t f = (BA0 * lambdaBA * (TMath::Exp(-lambdaBA * timeVar)));

	f += (CS0 * lambdaBA * lambdaCS * ((TMath::Exp(-lambdaCS * timeVar)) / ((lambdaBA-lambdaCS))));
	f += (CS0 * lambdaBA * lambdaCS * ((TMath::Exp(-lambdaBA * timeVar)) / ((lambdaCS-lambdaBA))));

	return f;
}

Double_t FitFunction::BADecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CS0 = par[0];
	Double_t lambdaCS = par[1];
	Double_t BA0 = par[2];
	Double_t lambdaBA = par[3];

	Double_t f = BA0 * (1.0 - TMath::Exp(-lambdaBA * timeVar));

	f += CS0 * lambdaBA * (1.0 - TMath::Exp(-lambdaCS * timeVar)) / (lambdaBA-lambdaCS));
	f += CS0 * lambdaCS * (1.0 - TMath::Exp(-lambdaBA * timeVar)) / (lambdaCS-lambdaBA));

	return f;
}

Double_t FitFunction::LADecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CS0 = par[0];
	Double_t lambdaCS = par[1];
	Double_t BA0 = par[2];
	Double_t lambdaBA = par[3];
	Double_t LA0 = par[4];
	Double_t lambdaLA = par[5];

	Double_t f = (LA0 * lambdaLA * (TMath::Exp(-lambdaLA * timeVar)));

	f += (CS0 * lambdaLA * lambdaCS * lambdaBA * ((TMath::Exp(-lambdaCS * timeVar)) / ((lambdaBA-lambdaCS)*(lambdaLA-lambdaCS))));
	f += (CS0 * lambdaLA * lambdaCS * lambdaBA * ((TMath::Exp(-lambdaBA * timeVar)) / ((lambdaCS-lambdaBA)*(lambdaLA-lambdaBA))));
	f += (CS0 * lambdaLA * lambdaCS * lambdaBA * ((TMath::Exp(-lambdaLA * timeVar)) / ((lambdaCS-lambdaLA)*(lambdaBA-lambdaLA))));

	f += (BA0 * lambdaLA * lambdaBA * ((TMath::Exp(-lambdaBA * timeVar)) / ((lambdaLA-lambdaBA))));
	f += (BA0 * lambdaLA * lambdaBA * ((TMath::Exp(-lambdaLA * timeVar)) / ((lambdaBA-lambdaLA))));

	return f;
}

Double_t FitFunction::LADecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CS0 = par[0];
	Double_t lambdaCS = par[1];
	Double_t BA0 = par[2];
	Double_t lambdaBA = par[3];
	Double_t LA0 = par[4];
	Double_t lambdaLA = par[5];

	Double_t f = LA0 * (1.0 - TMath::Exp(-lambdaLA * timeVar));

	f += CS0 * lambdaBA * lambdaLA * (1.0 - TMath::Exp(-lambdaCS * timeVar)) / (lambdaBA-lambdaCS)*(lambdaLA-lambdaCS));
	f += CS0 * lambdaCS * lambdaLA * (1.0 - TMath::Exp(-lambdaBA * timeVar)) / (lambdaCS-lambdaBA)*(lambdaLA-lambdaBA));
	f += CS0 * lambdaCS * lambdaBA * (1.0 - TMath::Exp(-lambdaLA * timeVar)) / (lambdaCS-lambdaLA)*(lambdaBA-lambdaLA));

	f += BA0 * lambdaLA * (1.0 - TMath::Exp(-lambdaBA * timeVar)) / (lambdaLA-lambdaBA));
	f += BA0 * lambdaBA * (1.0 - TMath::Exp(-lambdaLA * timeVar)) / (lambdaBA-lambdaLA));

	return f;
}

Double_t FitFunction::BatemanDecaybyActivty(Double_t *x, Double_t *par)
{
	Double_t hold = 0.0;
	for(int i = 0; i < numElements; i++)
	{
		hold += batemanFitFunctions[i](x, par);
	}
	return hold;
}

Double_t IntegralDecaybyActivity(Double_t *x, Double_t *par)
{
	Double_t hold = 0.0;
	for(int i = 0; i < numElements; i++)
	{
		hold += integralFitFunctions[i](x, par);
	}
	return hold;
}

endif;