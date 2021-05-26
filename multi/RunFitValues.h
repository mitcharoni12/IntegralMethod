/*
    PURPOSE: used to store arrays of base 4 run elements
    EX: fit of CS would have N0, N0Error, halfLife, halfLifeError. A 20x run would need arrays for store these 4 variables
    this class stores those
*/
#ifndef RUNFITVALUES_H
#define RUNFITVALUES_H

#include "GeneralValues.h"
#include "TMath.h"

class RunFitValues{
private:
    Int_t numRuns;
    GeneralValues* N0Run;
    GeneralValues* N0ErrorRun;
    GeneralValues* halfLifeRun;
    GeneralValues* halfLifeErrorRun;
public:
    RunFitValues(Int_t numRuns);
    ~RunFitValues();
    void SetAnN0(Int_t runIndex, Double_t N0);
    void SetAnN0Error(Int_t runIndex, Double_t N0Error);
    void SetAnHalfLife(Int_t runIndex, Double_t HalfLife);
    void SetAnHalfLifeError(Int_t runIndex, Double_t halfLifeError);
    Double_t GetAnN0(Int_t runIndex);
    Double_t GetAnN0Error(Int_t runIndex);
    Double_t GetAnHalfLife(Int_t runIndex);
    Double_t GetAnHalfLifeError(Int_t runIndex);
    Double_t* GetN0Arr(){return N0Run->GetFitArr();}
    Double_t* GetN0ErrorArr(){return N0ErrorRun->GetFitArr();}
    Double_t* GetHalfLifeArr(){return halfLifeRun->GetFitArr();}
    Double_t* GetHalfLifeErrorArr(){return halfLifeErrorRun->GetFitArr();}
};

RunFitValues::RunFitValues(Int_t numRuns)
{
    this->numRuns = numRuns;
    N0Run = new GeneralValues(numRuns);
    N0ErrorRun = new GeneralValues(numRuns);
    halfLifeRun = new GeneralValues(numRuns);
    halfLifeErrorRun = new GeneralValues(numRuns);
}

RunFitValues::~RunFitValues()
{
    delete N0Run;
    delete N0ErrorRun;
    delete halfLifeRun;
    delete halfLifeErrorRun;
}

void RunFitValues::SetAnN0(Int_t runIndex, Double_t N0)
{
    N0Run->SetVal(runIndex, N0);
}

void RunFitValues::SetAnN0Error(Int_t runIndex, Double_t N0Error)
{
    N0ErrorRun->SetVal(runIndex, N0Error);
}

void RunFitValues::SetAnHalfLife(Int_t runIndex, Double_t HalfLife)
{
    halfLifeRun->SetVal(runIndex, HalfLife);
}

void RunFitValues::SetAnHalfLifeError(Int_t runIndex, Double_t halfLifeError)
{
    halfLifeErrorRun->SetVal(runIndex, halfLifeError);
}

Double_t RunFitValues::GetAnN0(Int_t runIndex)
{
    return N0Run->GetVal(runIndex);
}

Double_t RunFitValues::GetAnN0Error(Int_t runIndex)
{
    return N0ErrorRun->GetVal(runIndex);
}

Double_t RunFitValues::GetAnHalfLife(Int_t runIndex)
{
    return halfLifeRun->GetVal(runIndex);
}

Double_t RunFitValues::GetAnHalfLifeError(Int_t runIndex)
{
    return halfLifeErrorRun->GetVal(runIndex);
}

#endif