/*
    PURPOSE: stores run values for a decay chain
    EX: each element in decay chain fit of La->Ba->Cs has elements N0, N0Error, halfLife, halfLifeError. This class stores an array
    of RunFitValues
*/
#ifndef CHAINRUNFITVALUES_H
#define CHAINRUNFITVALUES_H

#include "RunFitValues.h"
#include "TMath.h"

class ChainRunFitValues{
private:
    Int_t numChainElements;
    Int_t numRuns;
    RunFitValues** chainRunValues;
public:
    ChainRunFitValues(Int_t numChainElements, Int_t numRuns);
    ~ChainRunFitValues();
    void SetAnN0(Int_t runIndex, Int_t elementIndex, Double_t N0);
    void SetAnN0Error(Int_t runIndex, Int_t elementIndex, Double_t N0Error);
    void SetAnHalfLife(Int_t runIndex, Int_t elementIndex, Double_t halfLife);
    void SetAnHalfLifeError(Int_t runIndex, Int_t elementIndex, Double_t halfLifeError);
    Double_t GetAnN0(Int_t runIndex, Int_t elementIndex);
    Double_t GetAnN0Error(Int_t runIndex, Int_t elementIndex);
    Double_t GetAnHalfLife(Int_t runIndex, Int_t elementIndex);
    Double_t GetAnHalfLifeError(Int_t runIndex, Int_t elementIndex);
    Double_t* GetN0Arr(Int_t elementIndex);
    Double_t* GetN0ErrorArr(Int_t elementIndex);
    Double_t* GetHalfLifeArr(Int_t elementIndex);
    Double_t* GetHalfLifeErrorArr(Int_t elementIndex);
};

ChainRunFitValues::ChainRunFitValues(Int_t numChainElements, Int_t numRuns)
{
    this->numChainElements = numChainElements;
    this->numRuns = numRuns;
    chainRunValues = new RunFitValues* [numChainElements];
    for(int i = 0; i < numChainElements; i++)
    {
        chainRunValues[i] = new RunFitValues(numRuns);
    }
}

ChainRunFitValues::~ChainRunFitValues()
{
    for(int i = 0; i < numChainElements; i++)
    {
        delete chainRunValues[i];
    }
    delete [] chainRunValues;
}

void ChainRunFitValues::SetAnN0(Int_t runIndex, Int_t elementIndex, Double_t N0)
{
    chainRunValues[elementIndex]->SetAnN0(runIndex, N0);
}

void ChainRunFitValues::SetAnN0Error(Int_t runIndex, Int_t elementIndex, Double_t N0Error)
{
    chainRunValues[elementIndex]->SetAnN0Error(runIndex, N0Error);
}

void ChainRunFitValues::SetAnHalfLife(Int_t runIndex, Int_t elementIndex, Double_t halfLife)
{
    chainRunValues[elementIndex]->SetAnHalfLife(runIndex, halfLife);
}

void ChainRunFitValues::SetAnHalfLifeError(Int_t runIndex, Int_t elementIndex, Double_t halfLifeError)
{
    chainRunValues[elementIndex]->SetAnHalfLifeError(runIndex, halfLifeError);
}

Double_t ChainRunFitValues::GetAnN0(Int_t runIndex, Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetAnN0(runIndex);
}

Double_t ChainRunFitValues::GetAnN0Error(Int_t runIndex, Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetAnN0Error(runIndex);
}

Double_t ChainRunFitValues::GetAnHalfLife(Int_t runIndex, Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetAnHalfLife(runIndex);
}

Double_t ChainRunFitValues::GetAnHalfLifeError(Int_t runIndex, Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetAnHalfLifeError(runIndex);
}

Double_t* ChainRunFitValues::GetN0Arr(Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetN0Arr();
}

Double_t* ChainRunFitValues::GetN0ErrorArr(Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetN0ErrorArr();
}

Double_t* ChainRunFitValues::GetHalfLifeArr(Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetHalfLifeArr();
}

Double_t* ChainRunFitValues::GetHalfLifeErrorArr(Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetHalfLifeErrorArr();
}

#endif