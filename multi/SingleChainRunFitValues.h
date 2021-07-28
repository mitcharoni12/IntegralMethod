#ifndef SINGLECHAINRUNFITVALUES_H
#define SINGLECHAINRUNFITVALUES_H

#include "ChainRunFitValues.h"
#include "TMath.h"

/// Used to store all the fit values for each element in the Bateman/integral histograms
/// EX: in the chain 144Cs->144Ba->144La there are single Bateman histograms for deacys of Cs, Ba and La. However, the Bateman equation for Ba contains componets for Ba and Cs and the Bateman equation or La contains componsets for Cs, Ba, and La.
/// Each componet in each element of the chain each has values fitted values for N0, N0 error, half life, and half life error.
/// This class stores the fitted values for Cs which just has fitted values for Cs, the fitted values Ba which has fitted values for Cs and Ba, and the fitted values for La which has Cs, Ba and La.
class SingleChainRunFitValues{
private:
    Int_t numChainElements;
    Int_t numRuns;
    ChainRunFitValues** singleChainRunVals;
public:
    SingleChainRunFitValues(Int_t numChainElements, Int_t numRuns);
    ~SingleChainRunFitValues();
    void SetAnN0(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex, Double_t N0);
    void SetAnN0Error(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex, Double_t N0Error);
    void SetAnHalfLife(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex, Double_t halfLife);
    void SetAnHalfLifeError(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex, Double_t halfLifeError);
    Double_t GetAnN0(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex);
    Double_t GetAnN0Error(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex);
    Double_t GetAnHalfLife(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex);
    Double_t GetAnHalfLifeError(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex);
    Double_t* GetN0Arr(Int_t elementIndex, Int_t elementSubIndex);
    Double_t* GetN0ErrorArr(Int_t elementIndex, Int_t elementSubIndex);
    Double_t* GetHalfLifeArr(Int_t elementIndex, Int_t elementSubIndex);
    Double_t* GetHalfLifeErrorArr(Int_t elementIndex, Int_t elementSubIndex);
};

SingleChainRunFitValues::SingleChainRunFitValues(Int_t numChainElements, Int_t numRuns)
{
    this->numChainElements = numChainElements;
    this->numRuns = numRuns;
    singleChainRunVals = new ChainRunFitValues* [numChainElements];
    for(int i = 0; i < numChainElements; i++)
    {
        singleChainRunVals[i] = new ChainRunFitValues(i+1, numRuns);
    }
}

SingleChainRunFitValues::~SingleChainRunFitValues()
{
    for(int i = 0; i < numChainElements; i++)
    {
        delete singleChainRunVals[i];
    }
    delete [] singleChainRunVals;
}

/// Sets an N0 at a specific run, element and element sub index.
void SingleChainRunFitValues::SetAnN0(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex, Double_t N0)
{
    singleChainRunVals[elementIndex]->SetAnN0(runIndex, elementSubIndex, N0);
}

/// Sets an N0Error at a specific run, element and element sub index.
void SingleChainRunFitValues::SetAnN0Error(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex, Double_t N0Error)
{
    singleChainRunVals[elementIndex]->SetAnN0Error(runIndex, elementSubIndex, N0Error);
}

/// Sets an Half life at a specific run, element and element sub index.
void SingleChainRunFitValues::SetAnHalfLife(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex, Double_t halfLife)
{
    singleChainRunVals[elementIndex]->SetAnHalfLife(runIndex, elementSubIndex, halfLife);
}

/// Sets an Half life error at a specific run, element and element sub index.
void SingleChainRunFitValues::SetAnHalfLifeError(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex, Double_t halfLifeError)
{
    singleChainRunVals[elementIndex]->SetAnHalfLifeError(runIndex, elementSubIndex, halfLifeError);
}

/// Gets an N0 at a specific run, element and element sub index.
Double_t SingleChainRunFitValues::GetAnN0(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex)
{
    return singleChainRunVals[elementIndex]->GetAnN0(runIndex, elementSubIndex);
}

/// Gets an N0Error at a specific run, element and element sub index.
Double_t SingleChainRunFitValues::GetAnN0Error(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex)
{
    return singleChainRunVals[elementIndex]->GetAnN0Error(runIndex, elementSubIndex);
}

/// Gets an Half life at a specific run, element and element sub index.
Double_t SingleChainRunFitValues::GetAnHalfLife(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex)
{
    return singleChainRunVals[elementIndex]->GetAnHalfLife(runIndex, elementSubIndex);
}

/// Gets an Half life error at a specific run, element and element sub index.
Double_t SingleChainRunFitValues::GetAnHalfLifeError(Int_t runIndex, Int_t elementIndex, Int_t elementSubIndex)
{
    return singleChainRunVals[elementIndex]->GetAnHalfLifeError(runIndex, elementSubIndex);
}

/// Gets the array of N0 values for a cycle at a specific element and element sub index.
Double_t* SingleChainRunFitValues::GetN0Arr(Int_t elementIndex, Int_t elementSubIndex)
{
    return singleChainRunVals[elementIndex]->GetN0Arr(elementSubIndex);
}

/// Gets the array of N0 error values for a cycle at a specific element and element sub index.
Double_t* SingleChainRunFitValues::GetN0ErrorArr(Int_t elementIndex, Int_t elementSubIndex)
{
    return singleChainRunVals[elementIndex]->GetN0ErrorArr(elementSubIndex);
}

/// Gets the array of Half life values for a cycle at a specific element and element sub index.
Double_t* SingleChainRunFitValues::GetHalfLifeArr(Int_t elementIndex, Int_t elementSubIndex)
{
    return singleChainRunVals[elementIndex]->GetHalfLifeArr(elementSubIndex);
}

/// Gets the array of Half life error values for a cycle at a specific element and element sub index.
Double_t* SingleChainRunFitValues::GetHalfLifeErrorArr(Int_t elementIndex, Int_t elementSubIndex)
{
    return singleChainRunVals[elementIndex]->GetHalfLifeErrorArr(elementSubIndex);
}

#endif