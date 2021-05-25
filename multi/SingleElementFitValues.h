/*
    PURPOSE: Way to store fitted values for each element in decay chain individually
        EX: in decay chain 144CS->144Ba->144LA. The formula for the decay of JUST LA has componets of CS, BA, and LA. we need a way to store single fit values for
        EACH ELEMENT IN THE DECAY CHAIN
*/
#ifndef SINGLEELEMENTFITVALUES_H
#define SINGLEELEMENTFITVALUES_H

#include "ChainFitValues.h"

class SingleElementFitValues{
private:
    ChainFitValues** singleFitValues;
    Int_t numChainElements;
public:
    void SetElementChainFitValues(ChainFitValues** singleFitValues){this->singleFitValues = singleFitValues;}
    ChainFitValues** GetElementChainFitValues(){return singleFitValues;}
    ChainFitValues* GetChainFitValues(Int_t i){return singleFitValues[i];}
    void SetAnN0(int chainNum, int elementNum, Double_t N0);
    void SetAnN0Error(int chainNum, int elementNum, Double_t N0Error);
    void SetAnHalfLife(int chainNum, int elementNum, Double_t halfLife);
    void SetAnHalfLifeError(int chainNum, int elementNum, Double_t halfLifeError);
    Double_t GetAnN0(int chainNum, int elementNum);
    Double_t GetAnN0Error(int chainNum, int elementNum);
    Double_t GetAnHalfLife(int chainNum, int elementNum);
    Double_t GetAnHalfLifeError(int chainNum, int elementNum);
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

void SingleElementFitValues::SetAnN0(int chainNum, int elementNum, Double_t N0)
{
    singleFitValues[chainNum]->SetAnN0(elementNum, N0);
}

void SingleElementFitValues::SetAnN0Error(int chainNum, int elementNum, Double_t N0Error)
{
    singleFitValues[chainNum]->SetAnN0Error(elementNum, N0Error);
}

void SingleElementFitValues::SetAnHalfLife(int chainNum, int elementNum, Double_t halfLife)
{  
    singleFitValues[chainNum]->SetAnHalfLife(elementNum, halfLife);
}

void SingleElementFitValues::SetAnHalfLifeError(int chainNum, int elementNum, Double_t halfLifeError)
{
    singleFitValues[chainNum]->SetAnHalfLifeError(elementNum, halfLifeError);
}

Double_t SingleElementFitValues::GetAnN0(int chainNum, int elementNum)
{
    Double_t N0 = singleFitValues[chainNum]->GetAnN0(elementNum);
    return N0;
}

Double_t SingleElementFitValues::GetAnN0Error(int chainNum, int elementNum)
{
    Double_t N0Error = singleFitValues[chainNum]->GetAnN0Error(elementNum);
    return N0Error;
}

Double_t SingleElementFitValues::GetAnHalfLife(int chainNum, int elementNum)
{
    Double_t halfLife = singleFitValues[chainNum]->GetAnHalfLife(elementNum);
    return halfLife;
}

Double_t SingleElementFitValues::GetAnHalfLifeError(int chainNum, int elementNum)
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
    delete singleFitValues;
}
#endif