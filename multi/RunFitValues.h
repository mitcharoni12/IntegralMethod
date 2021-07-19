#ifndef RUNFITVALUES_H
#define RUNFITVALUES_H

#include "GeneralValues.h"
#include "TMath.h"

/// Used for storing the the fitted values for either the values of a total Bateman/integral histogram or one element of a single Bateman/integral histogram of multiple runs.
/// EX: If we look at decay chain 144Cs->144Ba->144La and we look at the total Bateman histogram, the total bateman equation is going to have values for N0, N0Error, half life, and half life error for each element in that chain and in a situation with 20x runs we need to store these in an array type format.
// This class stores that type of data.
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

/// Sets an N0 for a specified run index.
void RunFitValues::SetAnN0(Int_t runIndex, Double_t N0)
{
    N0Run->SetVal(runIndex, N0);
}

/// Sets an N0 for a specified run index.
void RunFitValues::SetAnN0Error(Int_t runIndex, Double_t N0Error)
{
    N0ErrorRun->SetVal(runIndex, N0Error);
}

/// Sets an N0 for a specified run index.
void RunFitValues::SetAnHalfLife(Int_t runIndex, Double_t HalfLife)
{
    halfLifeRun->SetVal(runIndex, HalfLife);
}

/// Sets an N0 for a specified run index.
void RunFitValues::SetAnHalfLifeError(Int_t runIndex, Double_t halfLifeError)
{
    halfLifeErrorRun->SetVal(runIndex, halfLifeError);
}

/// Gets an N0 for a specified run index.
Double_t RunFitValues::GetAnN0(Int_t runIndex)
{
    return N0Run->GetVal(runIndex);
}

/// Gets an N0 error for a specified run index.
Double_t RunFitValues::GetAnN0Error(Int_t runIndex)
{
    return N0ErrorRun->GetVal(runIndex);
}

/// Gets an half life for a specified run index.
Double_t RunFitValues::GetAnHalfLife(Int_t runIndex)
{
    return halfLifeRun->GetVal(runIndex);
}

/// Gets an half life error for a specified run index.
Double_t RunFitValues::GetAnHalfLifeError(Int_t runIndex)
{
    return halfLifeErrorRun->GetVal(runIndex);
}

#endif