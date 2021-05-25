/*
    PURPOSE: way to store fitted values for a decay chain fit.
        EX: decay chain fit of CS144->BA144->LA144 has N0, N0Error, HalfLife, HalfLifeError for each element ion decay chain
*/
#ifndef CHAINFITVALUES_H
#define CHAINFITVALUES_H

#include "FitValues.h"

class ChainFitValues{
private:
    FitValues** chainFitVals;
    Int_t numElements;
public:
    void SetChainFitValues(FitValues** chainFitVals){this->chainFitVals = chainFitVals;}
    Int_t GetNumElement(){return numElements;}
    FitValues** GetChainFitValues(){return chainFitVals;}
    FitValues* GetSingleElementFitValues(Int_t i);
    void SetAnN0(int elementNum, Double_t N0);
    void SetAnN0Error(int elementNum, Double_t N0Error);
    void SetAnHalfLife(int elementNum, Double_t halfLife);
    void SetAnHalfLifeError(int elementNum, Double_t halfLifeError);
    Double_t GetAnN0(int elementNum);
    Double_t GetAnN0Error(int elementNum);
    Double_t GetAnHalfLife(int elementNum);
    Double_t GetAnHalfLifeError(int elementNum);
    
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

void ChainFitValues::SetAnN0(int elementNum, Double_t N0)
{
    chainFitVals[elementNum]->SetN0(N0);
}

void ChainFitValues::SetAnN0Error(int elementNum, Double_t N0Error)
{
    chainFitVals[elementNum]->SetN0Error(N0Error);
}

void ChainFitValues::SetAnHalfLife(int elementNum, Double_t halfLife)
{
    chainFitVals[elementNum]->SetHalfLife(halfLife);
}

void ChainFitValues::SetAnHalfLifeError(int elementNum, Double_t halfLifeError)
{
    chainFitVals[elementNum]->SetHalfLifeError(halfLifeError);
}

Double_t ChainFitValues::GetAnN0(int elementNum)
{
    Double_t N0 = chainFitVals[elementNum]->GetN0();
    return N0;
}

Double_t ChainFitValues::GetAnN0Error(int elementNum)
{
    Double_t N0Error = chainFitVals[elementNum]->GetN0Error();
    return N0Error;
}

Double_t ChainFitValues::GetAnHalfLife(int elementNum)
{
    Double_t halfLife = chainFitVals[elementNum]->GetHalfLife();
    return halfLife;
}

Double_t ChainFitValues::GetAnHalfLifeError(int elementNum)
{
    Double_t halfLifeError = chainFitVals[elementNum]->GetHalfLifeError();
    return halfLifeError;
}

//returns certain element fit values in chain
FitValues* ChainFitValues::GetSingleElementFitValues(Int_t i)
{
    if(i < numElements && i > -1)
    {
        return chainFitVals[i];
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
    delete chainFitVals;
}

#endif