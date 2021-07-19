#ifndef CHAINRUNFITVALUES_H
#define CHAINRUNFITVALUES_H

#include "RunFitValues.h"
#include "TMath.h"

/// Used for storing the the fitted values for either the values of a total Bateman/integral histogram or one element of a single Bateman/integral histogram for multiple runs.
/// EX: If we look at decay chain 144Cs->144Ba->144La and we look at the total Bateman histogram, the total bateman equation is going to have values for N0, N0Error, half life, and half life error for each element in that chain and in a 20x run it would have 20 of those.
// This class stores that type of data.
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

/// Sets an N0 for a specified run index and specified element index.
void ChainRunFitValues::SetAnN0(Int_t runIndex, Int_t elementIndex, Double_t N0)
{
    chainRunValues[elementIndex]->SetAnN0(runIndex, N0);
}

/// Sets an N0 error for a specified run index and specified element index.
void ChainRunFitValues::SetAnN0Error(Int_t runIndex, Int_t elementIndex, Double_t N0Error)
{
    chainRunValues[elementIndex]->SetAnN0Error(runIndex, N0Error);
}

/// Sets an half life for a specified run index and specified element index.
void ChainRunFitValues::SetAnHalfLife(Int_t runIndex, Int_t elementIndex, Double_t halfLife)
{
    chainRunValues[elementIndex]->SetAnHalfLife(runIndex, halfLife);
}

/// Sets an half life error for a specified run index and specified element index.
void ChainRunFitValues::SetAnHalfLifeError(Int_t runIndex, Int_t elementIndex, Double_t halfLifeError)
{
    chainRunValues[elementIndex]->SetAnHalfLifeError(runIndex, halfLifeError);
}

/// Gets an N0 for a specified run index and specified element index.
Double_t ChainRunFitValues::GetAnN0(Int_t runIndex, Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetAnN0(runIndex);
}

/// Gets an N0 error for a specified run index and specified element index.
Double_t ChainRunFitValues::GetAnN0Error(Int_t runIndex, Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetAnN0Error(runIndex);
}

/// Gets an half life for a specified run index and specified element index.
Double_t ChainRunFitValues::GetAnHalfLife(Int_t runIndex, Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetAnHalfLife(runIndex);
}

/// Gets an half life error for a specified run index and specified element index.
Double_t ChainRunFitValues::GetAnHalfLifeError(Int_t runIndex, Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetAnHalfLifeError(runIndex);
}

/// Gets the array of N0 for the multiple runs of a specified element.
Double_t* ChainRunFitValues::GetN0Arr(Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetN0Arr();
}

/// Gets the array of N0 error for the multiple runs of a specified element.
Double_t* ChainRunFitValues::GetN0ErrorArr(Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetN0ErrorArr();
}

/// Gets the array of half life for the multiple runs of a specified element.
Double_t* ChainRunFitValues::GetHalfLifeArr(Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetHalfLifeArr();
}

/// Gets the array of half life error for the multiple runs of a specified element.
Double_t* ChainRunFitValues::GetHalfLifeErrorArr(Int_t elementIndex)
{
    return chainRunValues[elementIndex]->GetHalfLifeErrorArr();
}

#endif