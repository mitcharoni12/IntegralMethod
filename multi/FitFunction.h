#ifndef FITFUNCTION_H
#define FITFUNCTION_H

#include "TMath.h"

using namespace std;

class FitFunction{
public:
	typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);
	static Double_t MGDecayByActivity(Double_t *x, Double_t* par);
	static Double_t MGDecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t ALDecayByActivity(Double_t *x, Double_t* par);
	static Double_t ALDecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t SIDecayByActivity(Double_t *x, Double_t* par);
	static Double_t SIDecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t PDecayByActivity(Double_t *x, Double_t* par);
	static Double_t PDecayByActivityIntegral(Double_t *x, Double_t* par);
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
	batemanFitFunctions[0] = MGDecayByActivity;
	integralFitFunctions[0] = MGDecayByActivityIntegral;
	batemanFitFunctions[1] = ALDecayByActivity;
	integralFitFunctions[1] = ALDecayByActivityIntegral;
	batemanFitFunctions[2] = SIDecayByActivity;
	integralFitFunctions[2] = SIDecayByActivityIntegral;
	batemanFitFunctions[3] = PDecayByActivity;
	integralFitFunctions[3] = PDecayByActivityIntegral;
}

Double_t FitFunction::MGDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t MG0 = par[0];
	Double_t lambdaMG = par[1];

	Double_t f = (MG0 * lambdaMG * (TMath::Exp(-lambdaMG * timeVar)));

	return f;
}

Double_t FitFunction::MGDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t MG0 = par[0];
	Double_t lambdaMG = par[1];

	Double_t f = MG0 * (1.0 - TMath::Exp(-lambdaMG * timeVar));

	return f;
}

Double_t FitFunction::ALDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t MG0 = par[0];
	Double_t lambdaMG = par[1];
	Double_t AL0 = par[2];
	Double_t lambdaAL = par[3];

	Double_t f = (AL0 * lambdaAL * (TMath::Exp(-lambdaAL * timeVar)));

	f += (MG0 * lambdaAL * lambdaMG * ((TMath::Exp(-lambdaMG * timeVar)) / ((lambdaAL-lambdaMG))));
	f += (MG0 * lambdaAL * lambdaMG * ((TMath::Exp(-lambdaAL * timeVar)) / ((lambdaMG-lambdaAL))));

	return f;
}

Double_t FitFunction::ALDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t MG0 = par[0];
	Double_t lambdaMG = par[1];
	Double_t AL0 = par[2];
	Double_t lambdaAL = par[3];

	Double_t f = AL0 * (1.0 - TMath::Exp(-lambdaAL * timeVar));

	f += MG0 * lambdaAL * (1.0 - TMath::Exp(-lambdaMG * timeVar)) / ((lambdaAL-lambdaMG));
	f += MG0 * lambdaMG * (1.0 - TMath::Exp(-lambdaAL * timeVar)) / ((lambdaMG-lambdaAL));

	return f;
}

Double_t FitFunction::SIDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t MG0 = par[0];
	Double_t lambdaMG = par[1];
	Double_t AL0 = par[2];
	Double_t lambdaAL = par[3];
	Double_t SI0 = par[4];
	Double_t lambdaSI = par[5];

	Double_t f = (SI0 * lambdaSI * (TMath::Exp(-lambdaSI * timeVar)));

	f += (MG0 * lambdaSI * lambdaMG * lambdaAL * ((TMath::Exp(-lambdaMG * timeVar)) / ((lambdaAL-lambdaMG)*(lambdaSI-lambdaMG))));
	f += (MG0 * lambdaSI * lambdaMG * lambdaAL * ((TMath::Exp(-lambdaAL * timeVar)) / ((lambdaMG-lambdaAL)*(lambdaSI-lambdaAL))));
	f += (MG0 * lambdaSI * lambdaMG * lambdaAL * ((TMath::Exp(-lambdaSI * timeVar)) / ((lambdaMG-lambdaSI)*(lambdaAL-lambdaSI))));

	f += (AL0 * lambdaSI * lambdaAL * ((TMath::Exp(-lambdaAL * timeVar)) / ((lambdaSI-lambdaAL))));
	f += (AL0 * lambdaSI * lambdaAL * ((TMath::Exp(-lambdaSI * timeVar)) / ((lambdaAL-lambdaSI))));

	return f;
}

Double_t FitFunction::SIDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t MG0 = par[0];
	Double_t lambdaMG = par[1];
	Double_t AL0 = par[2];
	Double_t lambdaAL = par[3];
	Double_t SI0 = par[4];
	Double_t lambdaSI = par[5];

	Double_t f = SI0 * (1.0 - TMath::Exp(-lambdaSI * timeVar));

	f += MG0 * lambdaAL * lambdaSI * (1.0 - TMath::Exp(-lambdaMG * timeVar)) / ((lambdaAL-lambdaMG)*(lambdaSI-lambdaMG));
	f += MG0 * lambdaMG * lambdaSI * (1.0 - TMath::Exp(-lambdaAL * timeVar)) / ((lambdaMG-lambdaAL)*(lambdaSI-lambdaAL));
	f += MG0 * lambdaMG * lambdaAL * (1.0 - TMath::Exp(-lambdaSI * timeVar)) / ((lambdaMG-lambdaSI)*(lambdaAL-lambdaSI));

	f += AL0 * lambdaSI * (1.0 - TMath::Exp(-lambdaAL * timeVar)) / ((lambdaSI-lambdaAL));
	f += AL0 * lambdaAL * (1.0 - TMath::Exp(-lambdaSI * timeVar)) / ((lambdaAL-lambdaSI));

	return f;
}

Double_t FitFunction::PDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t MG0 = par[0];
	Double_t lambdaMG = par[1];
	Double_t AL0 = par[2];
	Double_t lambdaAL = par[3];
	Double_t SI0 = par[4];
	Double_t lambdaSI = par[5];
	Double_t P0 = par[6];
	Double_t lambdaP = par[7];

	Double_t f = (P0 * lambdaP * (TMath::Exp(-lambdaP * timeVar)));

	f += (MG0 * lambdaP * lambdaMG * lambdaAL * lambdaSI * ((TMath::Exp(-lambdaMG * timeVar)) / ((lambdaAL-lambdaMG)*(lambdaSI-lambdaMG)*(lambdaP-lambdaMG))));
	f += (MG0 * lambdaP * lambdaMG * lambdaAL * lambdaSI * ((TMath::Exp(-lambdaAL * timeVar)) / ((lambdaMG-lambdaAL)*(lambdaSI-lambdaAL)*(lambdaP-lambdaAL))));
	f += (MG0 * lambdaP * lambdaMG * lambdaAL * lambdaSI * ((TMath::Exp(-lambdaSI * timeVar)) / ((lambdaMG-lambdaSI)*(lambdaAL-lambdaSI)*(lambdaP-lambdaSI))));
	f += (MG0 * lambdaP * lambdaMG * lambdaAL * lambdaSI * ((TMath::Exp(-lambdaP * timeVar)) / ((lambdaMG-lambdaP)*(lambdaAL-lambdaP)*(lambdaSI-lambdaP))));

	f += (AL0 * lambdaP * lambdaAL * lambdaSI * ((TMath::Exp(-lambdaAL * timeVar)) / ((lambdaSI-lambdaAL)*(lambdaP-lambdaAL))));
	f += (AL0 * lambdaP * lambdaAL * lambdaSI * ((TMath::Exp(-lambdaSI * timeVar)) / ((lambdaAL-lambdaSI)*(lambdaP-lambdaSI))));
	f += (AL0 * lambdaP * lambdaAL * lambdaSI * ((TMath::Exp(-lambdaP * timeVar)) / ((lambdaAL-lambdaP)*(lambdaSI-lambdaP))));

	f += (SI0 * lambdaP * lambdaSI * ((TMath::Exp(-lambdaSI * timeVar)) / ((lambdaP-lambdaSI))));
	f += (SI0 * lambdaP * lambdaSI * ((TMath::Exp(-lambdaP * timeVar)) / ((lambdaSI-lambdaP))));

	return f;
}

Double_t FitFunction::PDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t MG0 = par[0];
	Double_t lambdaMG = par[1];
	Double_t AL0 = par[2];
	Double_t lambdaAL = par[3];
	Double_t SI0 = par[4];
	Double_t lambdaSI = par[5];
	Double_t P0 = par[6];
	Double_t lambdaP = par[7];

	Double_t f = P0 * (1.0 - TMath::Exp(-lambdaP * timeVar));

	f += MG0 * lambdaAL * lambdaSI * lambdaP * (1.0 - TMath::Exp(-lambdaMG * timeVar)) / ((lambdaAL-lambdaMG)*(lambdaSI-lambdaMG)*(lambdaP-lambdaMG));
	f += MG0 * lambdaMG * lambdaSI * lambdaP * (1.0 - TMath::Exp(-lambdaAL * timeVar)) / ((lambdaMG-lambdaAL)*(lambdaSI-lambdaAL)*(lambdaP-lambdaAL));
	f += MG0 * lambdaMG * lambdaAL * lambdaP * (1.0 - TMath::Exp(-lambdaSI * timeVar)) / ((lambdaMG-lambdaSI)*(lambdaAL-lambdaSI)*(lambdaP-lambdaSI));
	f += MG0 * lambdaMG * lambdaAL * lambdaSI * (1.0 - TMath::Exp(-lambdaP * timeVar)) / ((lambdaMG-lambdaP)*(lambdaAL-lambdaP)*(lambdaSI-lambdaP));

	f += AL0 * lambdaSI * lambdaP * (1.0 - TMath::Exp(-lambdaAL * timeVar)) / ((lambdaSI-lambdaAL)*(lambdaP-lambdaAL));
	f += AL0 * lambdaAL * lambdaP * (1.0 - TMath::Exp(-lambdaSI * timeVar)) / ((lambdaAL-lambdaSI)*(lambdaP-lambdaSI));
	f += AL0 * lambdaAL * lambdaSI * (1.0 - TMath::Exp(-lambdaP * timeVar)) / ((lambdaAL-lambdaP)*(lambdaSI-lambdaP));

	f += SI0 * lambdaP * (1.0 - TMath::Exp(-lambdaSI * timeVar)) / ((lambdaP-lambdaSI));
	f += SI0 * lambdaSI * (1.0 - TMath::Exp(-lambdaP * timeVar)) / ((lambdaSI-lambdaP));

	return f;
}

#endif