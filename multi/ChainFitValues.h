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
    
    ChainFitValues(Int_t numElements);
    ~ChainFitValues();
};

ChainFitValues::ChainFitValues(Int_t numElements)
{
    this->numElements = numElements;    
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