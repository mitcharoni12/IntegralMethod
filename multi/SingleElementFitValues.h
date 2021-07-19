#ifndef SINGLEELEMENTFITVALUES_H
#define SINGLEELEMENTFITVALUES_H

#include "ChainFitValues.h"
#include "TH1D.h"

/// Used to store all the fit values for each element in the Bateman/integral histograms
/// EX: in the chain 144Cs->144Ba->144La there are single Bateman histograms for deacys of Cs, Ba and La. However, the Bateman equation for Ba contains componets for Ba and Cs and the Bateman equation or La contains componsets for Cs, Ba, and La.
/// Each componet in each element of the chain each has values fitted values for N0, N0 error, half life, and half life error.
/// This class stores the fitted values for Cs which just has fitted values for Cs, the fitted values Ba which has fitted values for Cs and Ba, and the fitted values for La which has Cs, Ba and La.
class SingleElementFitValues{
private:
    ChainFitValues** singleFitValues;
    Int_t numChainElements;
public:
    void SetElementChainFitValues(ChainFitValues** singleFitValues){this->singleFitValues = singleFitValues;}
    ChainFitValues** GetElementChainFitValues(){return singleFitValues;}
    ChainFitValues* GetChainFitValues(Int_t i){return singleFitValues[i];}
    void SetAnN0(Int_t chainNum, Int_t elementNum, Double_t N0);
    void SetAnN0Error(Int_t chainNum, Int_t elementNum, Double_t N0Error);
    void SetAnHalfLife(Int_t chainNum, Int_t elementNum, Double_t halfLife);
    void SetAnHalfLifeError(Int_t chainNum, Int_t elementNum, Double_t halfLifeError);
    Double_t GetAnN0(Int_t chainNum, Int_t elementNum);
    Double_t GetAnN0Error(Int_t chainNum, Int_t elementNum);
    Double_t GetAnHalfLife(Int_t chainNum, Int_t elementNum);
    Double_t GetAnHalfLifeError(Int_t chainNum, Int_t elementNum);
    Int_t GetNumElements(){return numChainElements;}
    
    SingleElementFitValues(Int_t numChainElements);
    ~SingleElementFitValues();
};


SingleElementFitValues::SingleElementFitValues(Int_t numChainElements)
{
    this->numChainElements = numChainElements;
    singleFitValues = new ChainFitValues* [numChainElements];
    for(int i = 0; i < numChainElements; i++)
    {
        singleFitValues[i] = new ChainFitValues(i+1);
    }
}

/// Sets N0 for a specified componet of a specified element.
void SingleElementFitValues::SetAnN0(Int_t chainNum, Int_t elementNum, Double_t N0)
{
    singleFitValues[chainNum]->SetAnN0(elementNum, N0);
}

/// Sets N0 error for a specified componet of a specified element.
void SingleElementFitValues::SetAnN0Error(Int_t chainNum, Int_t elementNum, Double_t N0Error)
{
    singleFitValues[chainNum]->SetAnN0Error(elementNum, N0Error);
}

/// Sets half life for a specified componet of a specified element.
void SingleElementFitValues::SetAnHalfLife(Int_t chainNum, Int_t elementNum, Double_t halfLife)
{  
    singleFitValues[chainNum]->SetAnHalfLife(elementNum, halfLife);
}

/// Sets half life error for a specified componet of a specified element.
void SingleElementFitValues::SetAnHalfLifeError(Int_t chainNum, Int_t elementNum, Double_t halfLifeError)
{
    singleFitValues[chainNum]->SetAnHalfLifeError(elementNum, halfLifeError);
}

/// Gets N0 for a specified componet of a specified element.
Double_t SingleElementFitValues::GetAnN0(Int_t chainNum, Int_t elementNum)
{
    Double_t N0 = singleFitValues[chainNum]->GetAnN0(elementNum);
    return N0;
}

/// Gets N0 error for a specified componet of a specified element.
Double_t SingleElementFitValues::GetAnN0Error(Int_t chainNum, Int_t elementNum)
{
    Double_t N0Error = singleFitValues[chainNum]->GetAnN0Error(elementNum);
    return N0Error;
}

/// Gets half life for a specified componet of a specified element.
Double_t SingleElementFitValues::GetAnHalfLife(Int_t chainNum, Int_t elementNum)
{
    Double_t halfLife = singleFitValues[chainNum]->GetAnHalfLife(elementNum);
    return halfLife;
}

/// Gets half life for a specified componet of a specified element.
Double_t SingleElementFitValues::GetAnHalfLifeError(Int_t chainNum, Int_t elementNum)
{
    Double_t halfLifeError = singleFitValues[chainNum]->GetAnHalfLifeError(elementNum);
    return halfLifeError;
}

SingleElementFitValues::~SingleElementFitValues()
{
    for(int i = 0; i < numChainElements; i++)
    {
        delete singleFitValues[i];
    }
    delete [] singleFitValues;
}
#endif