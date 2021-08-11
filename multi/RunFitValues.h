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
    GeneralValues* N0Cycle;
    GeneralValues* N0ErrorCycle;
    GeneralValues* halfLifeCycle;
    GeneralValues* halfLifeErrorCycle;
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
    Double_t* GetN0Arr(){return N0Cycle->GetFitArr();}
    Double_t* GetN0ErrorArr(){return N0ErrorCycle->GetFitArr();}
    Double_t* GetHalfLifeArr(){return halfLifeCycle->GetFitArr();}
    Double_t* GetHalfLifeErrorArr(){return halfLifeErrorCycle->GetFitArr();}
};

RunFitValues::RunFitValues(Int_t numRuns)
{
    this->numRuns = numRuns;
    N0Cycle = new GeneralValues(numRuns);
    N0ErrorCycle = new GeneralValues(numRuns);
    halfLifeCycle = new GeneralValues(numRuns);
    halfLifeErrorCycle = new GeneralValues(numRuns);
}

RunFitValues::~RunFitValues()
{
    delete N0Cycle;
    delete N0ErrorCycle;
    delete halfLifeCycle;
    delete halfLifeErrorCycle;
}

/// Sets an N0 for a specified run index.
void RunFitValues::SetAnN0(Int_t runIndex, Double_t N0)
{
    N0Cycle->SetVal(runIndex, N0);
}

/// Sets an N0 for a specified run index.
void RunFitValues::SetAnN0Error(Int_t runIndex, Double_t N0Error)
{
    N0ErrorCycle->SetVal(runIndex, N0Error);
}

/// Sets an N0 for a specified run index.
void RunFitValues::SetAnHalfLife(Int_t runIndex, Double_t HalfLife)
{
    halfLifeCycle->SetVal(runIndex, HalfLife);
}

/// Sets an N0 for a specified run index.
void RunFitValues::SetAnHalfLifeError(Int_t runIndex, Double_t halfLifeError)
{
    halfLifeErrorCycle->SetVal(runIndex, halfLifeError);
}

/// Gets an N0 for a specified run index.
Double_t RunFitValues::GetAnN0(Int_t runIndex)
{
    return N0Cycle->GetVal(runIndex);
}

/// Gets an N0 error for a specified run index.
Double_t RunFitValues::GetAnN0Error(Int_t runIndex)
{
    return N0ErrorCycle->GetVal(runIndex);
}

/// Gets an half life for a specified run index.
Double_t RunFitValues::GetAnHalfLife(Int_t runIndex)
{
    return halfLifeCycle->GetVal(runIndex);
}

/// Gets an half life error for a specified run index.
Double_t RunFitValues::GetAnHalfLifeError(Int_t runIndex)
{
    return halfLifeErrorCycle->GetVal(runIndex);
}

#endif