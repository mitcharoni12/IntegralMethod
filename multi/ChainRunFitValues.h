#ifndef CHAINRUNFITVALUES_H
#define CHAINRUNFITVALUES_H

#include "GeneralValues.h"

class ChainRunFitValues{
private:
    GeneralValues** chainFitValues;
    int numChainElements
public:
    
    
    ChainRunFitValues(Int_t numChainElements)
    ~ChainRunFitValues();
};

ChainRunFitValues::ChainRunFitValues(Int_t numChainElements)
{
    this->numChainElements = numChainElements;
}

ChainRunFitValues::~ChainRunFitValues()
{
    for(int i = 0; i < numChainElements; i++)
    {
        delete chainFitValues[i];
    }
    delete chainFitValues;
}

#endif