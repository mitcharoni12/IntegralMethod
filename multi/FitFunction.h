#ifndef FITFUNCTION_H
#define FITFUNCTION_H

#include "TMath.h"

using namespace std;

class FitFunction{
public:
	typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);
	static Double_t CADecayByActivity(Double_t *x, Double_t* par);
	static Double_t CADecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t SCDecayByActivity(Double_t *x, Double_t* par);
	static Double_t SCDecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t TIDecayByActivity(Double_t *x, Double_t* par);
	static Double_t TIDecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t VDecayByActivity(Double_t *x, Double_t* par);
	static Double_t VDecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t CrDecayByActivity(Double_t *x, Double_t* par);
	static Double_t CrDecayByActivityIntegral(Double_t *x, Double_t* par);
	static Double_t MnDecayByActivity(Double_t *x, Double_t* par);
	static Double_t MnDecayByActivityIntegral(Double_t *x, Double_t* par);
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
	batemanFitFunctions[0] = CADecayByActivity;
	integralFitFunctions[0] = CADecayByActivityIntegral;
	batemanFitFunctions[1] = SCDecayByActivity;
	integralFitFunctions[1] = SCDecayByActivityIntegral;
	batemanFitFunctions[2] = TIDecayByActivity;
	integralFitFunctions[2] = TIDecayByActivityIntegral;
	batemanFitFunctions[3] = VDecayByActivity;
	integralFitFunctions[3] = VDecayByActivityIntegral;
	batemanFitFunctions[4] = CrDecayByActivity;
	integralFitFunctions[4] = CrDecayByActivityIntegral;
	batemanFitFunctions[5] = MnDecayByActivity;
	integralFitFunctions[5] = MnDecayByActivityIntegral;
}

Double_t FitFunction::CADecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];

	Double_t f = (CA0 * lambdaCA * (TMath::Exp(-lambdaCA * timeVar)));

	return f;
}

Double_t FitFunction::CADecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];

	Double_t f = CA0 * (1.0 - TMath::Exp(-lambdaCA * timeVar));

	return f;
}

Double_t FitFunction::SCDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];

	Double_t f = (SC0 * lambdaSC * (TMath::Exp(-lambdaSC * timeVar)));

	f += (CA0 * lambdaSC * lambdaCA * ((TMath::Exp(-lambdaCA * timeVar)) / ((lambdaSC-lambdaCA))));
	f += (CA0 * lambdaSC * lambdaCA * ((TMath::Exp(-lambdaSC * timeVar)) / ((lambdaCA-lambdaSC))));

	return f;
}

Double_t FitFunction::SCDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];

	Double_t f = SC0 * (1.0 - TMath::Exp(-lambdaSC * timeVar));

	f += CA0 * lambdaSC * (1.0 - TMath::Exp(-lambdaCA * timeVar)) / (lambdaSC-lambdaCA);
	f += CA0 * lambdaCA * (1.0 - TMath::Exp(-lambdaSC * timeVar)) / (lambdaCA-lambdaSC);

	return f;
}

Double_t FitFunction::TIDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];
	Double_t TI0 = par[4];
	Double_t lambdaTI = par[5];

	Double_t f = (TI0 * lambdaTI * (TMath::Exp(-lambdaTI * timeVar)));

	f += (CA0 * lambdaTI * lambdaCA * lambdaSC * ((TMath::Exp(-lambdaCA * timeVar)) / ((lambdaSC-lambdaCA)*(lambdaTI-lambdaCA))));
	f += (CA0 * lambdaTI * lambdaCA * lambdaSC * ((TMath::Exp(-lambdaSC * timeVar)) / ((lambdaCA-lambdaSC)*(lambdaTI-lambdaSC))));
	f += (CA0 * lambdaTI * lambdaCA * lambdaSC * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaCA-lambdaTI)*(lambdaSC-lambdaTI))));

	f += (SC0 * lambdaTI * lambdaSC * ((TMath::Exp(-lambdaSC * timeVar)) / ((lambdaTI-lambdaSC))));
	f += (SC0 * lambdaTI * lambdaSC * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaSC-lambdaTI))));

	return f;
}

Double_t FitFunction::TIDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];
	Double_t TI0 = par[4];
	Double_t lambdaTI = par[5];

	Double_t f = TI0 * (1.0 - TMath::Exp(-lambdaTI * timeVar));

	f += CA0 * lambdaSC * lambdaTI * (1.0 - TMath::Exp(-lambdaCA * timeVar)) / (lambdaSC-lambdaCA)*(lambdaTI-lambdaCA);
	f += CA0 * lambdaCA * lambdaTI * (1.0 - TMath::Exp(-lambdaSC * timeVar)) / (lambdaCA-lambdaSC)*(lambdaTI-lambdaSC);
	f += CA0 * lambdaCA * lambdaSC * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaCA-lambdaTI)*(lambdaSC-lambdaTI);

	f += SC0 * lambdaTI * (1.0 - TMath::Exp(-lambdaSC * timeVar)) / (lambdaTI-lambdaSC);
	f += SC0 * lambdaSC * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaSC-lambdaTI);

	return f;
}

Double_t FitFunction::VDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];
	Double_t TI0 = par[4];
	Double_t lambdaTI = par[5];
	Double_t V0 = par[6];
	Double_t lambdaV = par[7];

	Double_t f = (V0 * lambdaV * (TMath::Exp(-lambdaV * timeVar)));

	f += (CA0 * lambdaV * lambdaCA * lambdaSC * lambdaTI * ((TMath::Exp(-lambdaCA * timeVar)) / ((lambdaSC-lambdaCA)*(lambdaTI-lambdaCA)*(lambdaV-lambdaCA))));
	f += (CA0 * lambdaV * lambdaCA * lambdaSC * lambdaTI * ((TMath::Exp(-lambdaSC * timeVar)) / ((lambdaCA-lambdaSC)*(lambdaTI-lambdaSC)*(lambdaV-lambdaSC))));
	f += (CA0 * lambdaV * lambdaCA * lambdaSC * lambdaTI * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaCA-lambdaTI)*(lambdaSC-lambdaTI)*(lambdaV-lambdaTI))));
	f += (CA0 * lambdaV * lambdaCA * lambdaSC * lambdaTI * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaCA-lambdaV)*(lambdaSC-lambdaV)*(lambdaTI-lambdaV))));

	f += (SC0 * lambdaV * lambdaSC * lambdaTI * ((TMath::Exp(-lambdaSC * timeVar)) / ((lambdaTI-lambdaSC)*(lambdaV-lambdaSC))));
	f += (SC0 * lambdaV * lambdaSC * lambdaTI * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaSC-lambdaTI)*(lambdaV-lambdaTI))));
	f += (SC0 * lambdaV * lambdaSC * lambdaTI * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaSC-lambdaV)*(lambdaTI-lambdaV))));

	f += (TI0 * lambdaV * lambdaTI * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaV-lambdaTI))));
	f += (TI0 * lambdaV * lambdaTI * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaTI-lambdaV))));

	return f;
}

Double_t FitFunction::VDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];
	Double_t TI0 = par[4];
	Double_t lambdaTI = par[5];
	Double_t V0 = par[6];
	Double_t lambdaV = par[7];

	Double_t f = V0 * (1.0 - TMath::Exp(-lambdaV * timeVar));

	f += CA0 * lambdaSC * lambdaTI * lambdaV * (1.0 - TMath::Exp(-lambdaCA * timeVar)) / (lambdaSC-lambdaCA)*(lambdaTI-lambdaCA)*(lambdaV-lambdaCA);
	f += CA0 * lambdaCA * lambdaTI * lambdaV * (1.0 - TMath::Exp(-lambdaSC * timeVar)) / (lambdaCA-lambdaSC)*(lambdaTI-lambdaSC)*(lambdaV-lambdaSC);
	f += CA0 * lambdaCA * lambdaSC * lambdaV * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaCA-lambdaTI)*(lambdaSC-lambdaTI)*(lambdaV-lambdaTI);
	f += CA0 * lambdaCA * lambdaSC * lambdaTI * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaCA-lambdaV)*(lambdaSC-lambdaV)*(lambdaTI-lambdaV);

	f += SC0 * lambdaTI * lambdaV * (1.0 - TMath::Exp(-lambdaSC * timeVar)) / (lambdaTI-lambdaSC)*(lambdaV-lambdaSC);
	f += SC0 * lambdaSC * lambdaV * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaSC-lambdaTI)*(lambdaV-lambdaTI);
	f += SC0 * lambdaSC * lambdaTI * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaSC-lambdaV)*(lambdaTI-lambdaV);

	f += TI0 * lambdaV * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaV-lambdaTI);
	f += TI0 * lambdaTI * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaTI-lambdaV);

	return f;
}

Double_t FitFunction::CrDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];
	Double_t TI0 = par[4];
	Double_t lambdaTI = par[5];
	Double_t V0 = par[6];
	Double_t lambdaV = par[7];
	Double_t Cr0 = par[8];
	Double_t lambdaCr = par[9];

	Double_t f = (Cr0 * lambdaCr * (TMath::Exp(-lambdaCr * timeVar)));

	f += (CA0 * lambdaCr * lambdaCA * lambdaSC * lambdaTI * lambdaV * ((TMath::Exp(-lambdaCA * timeVar)) / ((lambdaSC-lambdaCA)*(lambdaTI-lambdaCA)*(lambdaV-lambdaCA)*(lambdaCr-lambdaCA))));
	f += (CA0 * lambdaCr * lambdaCA * lambdaSC * lambdaTI * lambdaV * ((TMath::Exp(-lambdaSC * timeVar)) / ((lambdaCA-lambdaSC)*(lambdaTI-lambdaSC)*(lambdaV-lambdaSC)*(lambdaCr-lambdaSC))));
	f += (CA0 * lambdaCr * lambdaCA * lambdaSC * lambdaTI * lambdaV * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaCA-lambdaTI)*(lambdaSC-lambdaTI)*(lambdaV-lambdaTI)*(lambdaCr-lambdaTI))));
	f += (CA0 * lambdaCr * lambdaCA * lambdaSC * lambdaTI * lambdaV * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaCA-lambdaV)*(lambdaSC-lambdaV)*(lambdaTI-lambdaV)*(lambdaCr-lambdaV))));
	f += (CA0 * lambdaCr * lambdaCA * lambdaSC * lambdaTI * lambdaV * ((TMath::Exp(-lambdaCr * timeVar)) / ((lambdaCA-lambdaCr)*(lambdaSC-lambdaCr)*(lambdaTI-lambdaCr)*(lambdaV-lambdaCr))));

	f += (SC0 * lambdaCr * lambdaSC * lambdaTI * lambdaV * ((TMath::Exp(-lambdaSC * timeVar)) / ((lambdaTI-lambdaSC)*(lambdaV-lambdaSC)*(lambdaCr-lambdaSC))));
	f += (SC0 * lambdaCr * lambdaSC * lambdaTI * lambdaV * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaSC-lambdaTI)*(lambdaV-lambdaTI)*(lambdaCr-lambdaTI))));
	f += (SC0 * lambdaCr * lambdaSC * lambdaTI * lambdaV * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaSC-lambdaV)*(lambdaTI-lambdaV)*(lambdaCr-lambdaV))));
	f += (SC0 * lambdaCr * lambdaSC * lambdaTI * lambdaV * ((TMath::Exp(-lambdaCr * timeVar)) / ((lambdaSC-lambdaCr)*(lambdaTI-lambdaCr)*(lambdaV-lambdaCr))));

	f += (TI0 * lambdaCr * lambdaTI * lambdaV * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaV-lambdaTI)*(lambdaCr-lambdaTI))));
	f += (TI0 * lambdaCr * lambdaTI * lambdaV * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaTI-lambdaV)*(lambdaCr-lambdaV))));
	f += (TI0 * lambdaCr * lambdaTI * lambdaV * ((TMath::Exp(-lambdaCr * timeVar)) / ((lambdaTI-lambdaCr)*(lambdaV-lambdaCr))));

	f += (V0 * lambdaCr * lambdaV * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaCr-lambdaV))));
	f += (V0 * lambdaCr * lambdaV * ((TMath::Exp(-lambdaCr * timeVar)) / ((lambdaV-lambdaCr))));

	return f;
}

Double_t FitFunction::CrDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];
	Double_t TI0 = par[4];
	Double_t lambdaTI = par[5];
	Double_t V0 = par[6];
	Double_t lambdaV = par[7];
	Double_t Cr0 = par[8];
	Double_t lambdaCr = par[9];

	Double_t f = Cr0 * (1.0 - TMath::Exp(-lambdaCr * timeVar));

	f += CA0 * lambdaSC * lambdaTI * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaCA * timeVar)) / (lambdaSC-lambdaCA)*(lambdaTI-lambdaCA)*(lambdaV-lambdaCA)*(lambdaCr-lambdaCA);
	f += CA0 * lambdaCA * lambdaTI * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaSC * timeVar)) / (lambdaCA-lambdaSC)*(lambdaTI-lambdaSC)*(lambdaV-lambdaSC)*(lambdaCr-lambdaSC);
	f += CA0 * lambdaCA * lambdaSC * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaCA-lambdaTI)*(lambdaSC-lambdaTI)*(lambdaV-lambdaTI)*(lambdaCr-lambdaTI);
	f += CA0 * lambdaCA * lambdaSC * lambdaTI * lambdaCr * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaCA-lambdaV)*(lambdaSC-lambdaV)*(lambdaTI-lambdaV)*(lambdaCr-lambdaV);
	f += CA0 * lambdaCA * lambdaSC * lambdaTI * lambdaV * (1.0 - TMath::Exp(-lambdaCr * timeVar)) / (lambdaCA-lambdaCr)*(lambdaSC-lambdaCr)*(lambdaTI-lambdaCr)*(lambdaV-lambdaCr);

	f += SC0 * lambdaTI * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaSC * timeVar)) / (lambdaTI-lambdaSC)*(lambdaV-lambdaSC)*(lambdaCr-lambdaSC);
	f += SC0 * lambdaSC * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaSC-lambdaTI)*(lambdaV-lambdaTI)*(lambdaCr-lambdaTI);
	f += SC0 * lambdaSC * lambdaTI * lambdaCr * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaSC-lambdaV)*(lambdaTI-lambdaV)*(lambdaCr-lambdaV);
	f += SC0 * lambdaSC * lambdaTI * lambdaV * (1.0 - TMath::Exp(-lambdaCr * timeVar)) / (lambdaSC-lambdaCr)*(lambdaTI-lambdaCr)*(lambdaV-lambdaCr);

	f += TI0 * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaV-lambdaTI)*(lambdaCr-lambdaTI);
	f += TI0 * lambdaTI * lambdaCr * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaTI-lambdaV)*(lambdaCr-lambdaV);
	f += TI0 * lambdaTI * lambdaV * (1.0 - TMath::Exp(-lambdaCr * timeVar)) / (lambdaTI-lambdaCr)*(lambdaV-lambdaCr);

	f += V0 * lambdaCr * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaCr-lambdaV);
	f += V0 * lambdaV * (1.0 - TMath::Exp(-lambdaCr * timeVar)) / (lambdaV-lambdaCr);

	return f;
}

Double_t FitFunction::MnDecayByActivity(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];
	Double_t TI0 = par[4];
	Double_t lambdaTI = par[5];
	Double_t V0 = par[6];
	Double_t lambdaV = par[7];
	Double_t Cr0 = par[8];
	Double_t lambdaCr = par[9];
	Double_t Mn0 = par[10];
	Double_t lambdaMn = par[11];

	Double_t f = (Mn0 * lambdaMn * (TMath::Exp(-lambdaMn * timeVar)));

	f += (CA0 * lambdaMn * lambdaCA * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaCA * timeVar)) / ((lambdaSC-lambdaCA)*(lambdaTI-lambdaCA)*(lambdaV-lambdaCA)*(lambdaCr-lambdaCA)*(lambdaMn-lambdaCA))));
	f += (CA0 * lambdaMn * lambdaCA * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaSC * timeVar)) / ((lambdaCA-lambdaSC)*(lambdaTI-lambdaSC)*(lambdaV-lambdaSC)*(lambdaCr-lambdaSC)*(lambdaMn-lambdaSC))));
	f += (CA0 * lambdaMn * lambdaCA * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaCA-lambdaTI)*(lambdaSC-lambdaTI)*(lambdaV-lambdaTI)*(lambdaCr-lambdaTI)*(lambdaMn-lambdaTI))));
	f += (CA0 * lambdaMn * lambdaCA * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaCA-lambdaV)*(lambdaSC-lambdaV)*(lambdaTI-lambdaV)*(lambdaCr-lambdaV)*(lambdaMn-lambdaV))));
	f += (CA0 * lambdaMn * lambdaCA * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaCr * timeVar)) / ((lambdaCA-lambdaCr)*(lambdaSC-lambdaCr)*(lambdaTI-lambdaCr)*(lambdaV-lambdaCr)*(lambdaMn-lambdaCr))));
	f += (CA0 * lambdaMn * lambdaCA * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaMn * timeVar)) / ((lambdaCA-lambdaMn)*(lambdaSC-lambdaMn)*(lambdaTI-lambdaMn)*(lambdaV-lambdaMn)*(lambdaCr-lambdaMn))));

	f += (SC0 * lambdaMn * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaSC * timeVar)) / ((lambdaTI-lambdaSC)*(lambdaV-lambdaSC)*(lambdaCr-lambdaSC)*(lambdaMn-lambdaSC))));
	f += (SC0 * lambdaMn * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaSC-lambdaTI)*(lambdaV-lambdaTI)*(lambdaCr-lambdaTI)*(lambdaMn-lambdaTI))));
	f += (SC0 * lambdaMn * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaSC-lambdaV)*(lambdaTI-lambdaV)*(lambdaCr-lambdaV)*(lambdaMn-lambdaV))));
	f += (SC0 * lambdaMn * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaCr * timeVar)) / ((lambdaSC-lambdaCr)*(lambdaTI-lambdaCr)*(lambdaV-lambdaCr)*(lambdaMn-lambdaCr))));
	f += (SC0 * lambdaMn * lambdaSC * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaMn * timeVar)) / ((lambdaSC-lambdaMn)*(lambdaTI-lambdaMn)*(lambdaV-lambdaMn)*(lambdaCr-lambdaMn))));

	f += (TI0 * lambdaMn * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaTI * timeVar)) / ((lambdaV-lambdaTI)*(lambdaCr-lambdaTI)*(lambdaMn-lambdaTI))));
	f += (TI0 * lambdaMn * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaTI-lambdaV)*(lambdaCr-lambdaV)*(lambdaMn-lambdaV))));
	f += (TI0 * lambdaMn * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaCr * timeVar)) / ((lambdaTI-lambdaCr)*(lambdaV-lambdaCr)*(lambdaMn-lambdaCr))));
	f += (TI0 * lambdaMn * lambdaTI * lambdaV * lambdaCr * ((TMath::Exp(-lambdaMn * timeVar)) / ((lambdaTI-lambdaMn)*(lambdaV-lambdaMn)*(lambdaCr-lambdaMn))));

	f += (V0 * lambdaMn * lambdaV * lambdaCr * ((TMath::Exp(-lambdaV * timeVar)) / ((lambdaCr-lambdaV)*(lambdaMn-lambdaV))));
	f += (V0 * lambdaMn * lambdaV * lambdaCr * ((TMath::Exp(-lambdaCr * timeVar)) / ((lambdaV-lambdaCr)*(lambdaMn-lambdaCr))));
	f += (V0 * lambdaMn * lambdaV * lambdaCr * ((TMath::Exp(-lambdaMn * timeVar)) / ((lambdaV-lambdaMn)*(lambdaCr-lambdaMn))));

	f += (Cr0 * lambdaMn * lambdaCr * ((TMath::Exp(-lambdaCr * timeVar)) / ((lambdaMn-lambdaCr))));
	f += (Cr0 * lambdaMn * lambdaCr * ((TMath::Exp(-lambdaMn * timeVar)) / ((lambdaCr-lambdaMn))));

	return f;
}

Double_t FitFunction::MnDecayByActivityIntegral(Double_t *x, Double_t *par)
{
	Float_t timeVar = x[0];
	Double_t CA0 = par[0];
	Double_t lambdaCA = par[1];
	Double_t SC0 = par[2];
	Double_t lambdaSC = par[3];
	Double_t TI0 = par[4];
	Double_t lambdaTI = par[5];
	Double_t V0 = par[6];
	Double_t lambdaV = par[7];
	Double_t Cr0 = par[8];
	Double_t lambdaCr = par[9];
	Double_t Mn0 = par[10];
	Double_t lambdaMn = par[11];

	Double_t f = Mn0 * (1.0 - TMath::Exp(-lambdaMn * timeVar));

	f += CA0 * lambdaSC * lambdaTI * lambdaV * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaCA * timeVar)) / (lambdaSC-lambdaCA)*(lambdaTI-lambdaCA)*(lambdaV-lambdaCA)*(lambdaCr-lambdaCA)*(lambdaMn-lambdaCA);
	f += CA0 * lambdaCA * lambdaTI * lambdaV * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaSC * timeVar)) / (lambdaCA-lambdaSC)*(lambdaTI-lambdaSC)*(lambdaV-lambdaSC)*(lambdaCr-lambdaSC)*(lambdaMn-lambdaSC);
	f += CA0 * lambdaCA * lambdaSC * lambdaV * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaCA-lambdaTI)*(lambdaSC-lambdaTI)*(lambdaV-lambdaTI)*(lambdaCr-lambdaTI)*(lambdaMn-lambdaTI);
	f += CA0 * lambdaCA * lambdaSC * lambdaTI * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaCA-lambdaV)*(lambdaSC-lambdaV)*(lambdaTI-lambdaV)*(lambdaCr-lambdaV)*(lambdaMn-lambdaV);
	f += CA0 * lambdaCA * lambdaSC * lambdaTI * lambdaV * lambdaMn * (1.0 - TMath::Exp(-lambdaCr * timeVar)) / (lambdaCA-lambdaCr)*(lambdaSC-lambdaCr)*(lambdaTI-lambdaCr)*(lambdaV-lambdaCr)*(lambdaMn-lambdaCr);
	f += CA0 * lambdaCA * lambdaSC * lambdaTI * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaMn * timeVar)) / (lambdaCA-lambdaMn)*(lambdaSC-lambdaMn)*(lambdaTI-lambdaMn)*(lambdaV-lambdaMn)*(lambdaCr-lambdaMn);

	f += SC0 * lambdaTI * lambdaV * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaSC * timeVar)) / (lambdaTI-lambdaSC)*(lambdaV-lambdaSC)*(lambdaCr-lambdaSC)*(lambdaMn-lambdaSC);
	f += SC0 * lambdaSC * lambdaV * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaSC-lambdaTI)*(lambdaV-lambdaTI)*(lambdaCr-lambdaTI)*(lambdaMn-lambdaTI);
	f += SC0 * lambdaSC * lambdaTI * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaSC-lambdaV)*(lambdaTI-lambdaV)*(lambdaCr-lambdaV)*(lambdaMn-lambdaV);
	f += SC0 * lambdaSC * lambdaTI * lambdaV * lambdaMn * (1.0 - TMath::Exp(-lambdaCr * timeVar)) / (lambdaSC-lambdaCr)*(lambdaTI-lambdaCr)*(lambdaV-lambdaCr)*(lambdaMn-lambdaCr);
	f += SC0 * lambdaSC * lambdaTI * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaMn * timeVar)) / (lambdaSC-lambdaMn)*(lambdaTI-lambdaMn)*(lambdaV-lambdaMn)*(lambdaCr-lambdaMn);

	f += TI0 * lambdaV * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaTI * timeVar)) / (lambdaV-lambdaTI)*(lambdaCr-lambdaTI)*(lambdaMn-lambdaTI);
	f += TI0 * lambdaTI * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaTI-lambdaV)*(lambdaCr-lambdaV)*(lambdaMn-lambdaV);
	f += TI0 * lambdaTI * lambdaV * lambdaMn * (1.0 - TMath::Exp(-lambdaCr * timeVar)) / (lambdaTI-lambdaCr)*(lambdaV-lambdaCr)*(lambdaMn-lambdaCr);
	f += TI0 * lambdaTI * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaMn * timeVar)) / (lambdaTI-lambdaMn)*(lambdaV-lambdaMn)*(lambdaCr-lambdaMn);

	f += V0 * lambdaCr * lambdaMn * (1.0 - TMath::Exp(-lambdaV * timeVar)) / (lambdaCr-lambdaV)*(lambdaMn-lambdaV);
	f += V0 * lambdaV * lambdaMn * (1.0 - TMath::Exp(-lambdaCr * timeVar)) / (lambdaV-lambdaCr)*(lambdaMn-lambdaCr);
	f += V0 * lambdaV * lambdaCr * (1.0 - TMath::Exp(-lambdaMn * timeVar)) / (lambdaV-lambdaMn)*(lambdaCr-lambdaMn);

	f += Cr0 * lambdaMn * (1.0 - TMath::Exp(-lambdaCr * timeVar)) / (lambdaMn-lambdaCr);
	f += Cr0 * lambdaCr * (1.0 - TMath::Exp(-lambdaMn * timeVar)) / (lambdaCr-lambdaMn);

	return f;
}

#endif