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
    Int_t GetNumElements(){return numChainElements;}
    
    SingleElementFitValues(Int_t numChainElements);
    ~SingleElementFitValues();
};


SingleElementFitValues::SingleElementFitValues(Int_t numChainElements)
{
    this->numChainElements = numChainElements;
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