#ifndef CHAINFITVALUES_H
#define CHAINFITVALUES_H

#include "FitValues.h"

using namespace std;

/// Used for storing the the fitted values for either the values of a total Bateman/integral histogram or one element of a single Bateman/integral histogram.
/// EX: If we look at decay chain 144Cs->144Ba->144La and we look at the total Bateman histogram, the total bateman equation is going to have values for N0, N0Error, half life, and half life error for each element in that chain.
/// This class stores that type of data.
class ChainFitValues{
private:
    FitValues** chainFitVals;
    Int_t numElements;
public:
    void SetChainFitValues(FitValues** chainFitVals){this->chainFitVals = chainFitVals;}
    Int_t GetNumElements(){return numElements;}
    FitValues** GetChainFitValues(){return chainFitVals;}
    FitValues* GetSingleElementFitValues(Int_t elementIndex);
    void SetAnN0(Int_t elementNum, Double_t N0);
    void SetAnN0Error(Int_t elementNum, Double_t N0Error);
    void SetAnHalfLife(Int_t elementNum, Double_t halfLife);
    void SetAnHalfLifeError(Int_t elementNum, Double_t halfLifeError);
    Double_t GetAnN0(Int_t elementNum);
    Double_t GetAnN0Error(Int_t elementNum);
    Double_t GetAnHalfLife(Int_t elementNum);
    Double_t GetAnHalfLifeError(Int_t elementNum);
    
    ChainFitValues(Int_t numElements);
    ~ChainFitValues();
};

ChainFitValues::ChainFitValues(Int_t numElements)
{
    this->numElements = numElements;
    chainFitVals = new FitValues* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        chainFitVals[i] = new FitValues();
    } 
}

/// Sets an N0 for a specified element
void ChainFitValues::SetAnN0(Int_t elementNum, Double_t N0)
{
    chainFitVals[elementNum]->SetN0(N0);
}

/// Sets an N0 error for a specified element
void ChainFitValues::SetAnN0Error(Int_t elementNum, Double_t N0Error)
{
    chainFitVals[elementNum]->SetN0Error(N0Error);
}

/// Sets an half life for a specified element
void ChainFitValues::SetAnHalfLife(Int_t elementNum, Double_t halfLife)
{
    chainFitVals[elementNum]->SetHalfLife(halfLife);
}

/// Sets an half life error for a specified element
void ChainFitValues::SetAnHalfLifeError(Int_t elementNum, Double_t halfLifeError)
{
    chainFitVals[elementNum]->SetHalfLifeError(halfLifeError);
}

/// Gets an N0 for a specified element
Double_t ChainFitValues::GetAnN0(Int_t elementNum)
{
    Double_t N0 = chainFitVals[elementNum]->GetN0();
    return N0;
}

/// Gets an N0 error for a specified element
Double_t ChainFitValues::GetAnN0Error(Int_t elementNum)
{
    Double_t N0Error = chainFitVals[elementNum]->GetN0Error();
    return N0Error;
}

/// Gets an half life for a specified element
Double_t ChainFitValues::GetAnHalfLife(Int_t elementNum)
{
    Double_t halfLife = chainFitVals[elementNum]->GetHalfLife();
    return halfLife;
}

/// Gets an half life error for a specified element
Double_t ChainFitValues::GetAnHalfLifeError(Int_t elementNum)
{
    Double_t halfLifeError = chainFitVals[elementNum]->GetHalfLifeError();
    return halfLifeError;
}

/// returns certain element fit values in chain
FitValues* ChainFitValues::GetSingleElementFitValues(Int_t elementIndex)
{
    if(elementIndex < numElements && elementIndex > -1)
    {
        return chainFitVals[elementIndex];
    }else{
        cout << "ERROR: invalid element index";
        return nullptr;
    }
}

ChainFitValues::~ChainFitValues()
{
    for(int i = 0; i < numElements; i++)
    {
        delete chainFitVals[i];
    }
    delete [] chainFitVals;
}

#endif