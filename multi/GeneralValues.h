/*
    PURPOSE: Basically just a glorified array, used to store array of values like N0 and half life for the entire durration of a run
    EX: 20x will need to store 20 different half life values
*/
#ifndef GENERALVALUES_H
#define GENERALVALUES_H

#include "TMath.h"

class GeneralValues{
private:
    Double_t* fitValues;
    Int_t numRuns;
public:
    GeneralValues(Int_t numRuns);
    ~GeneralValues();
    Int_t GetNumRuns(){return numRuns;}
    void SetVal(Int_t index, Double_t val);
    Double_t GetVal(Int_t index);
    Double_t* GetFitArr(){return fitValues;}
};

GeneralValues::GeneralValues(Int_t numRuns)
{
    this->numRuns = numRuns;
    fitValues = new Double_t[numRuns];
}

GeneralValues::~GeneralValues()
{
    delete [] fitValues;
}

//sets values in array fitValues
void GeneralValues::SetVal(Int_t index, Double_t val)
{
    fitValues[index] = val;
}

//gets values in array fitValues
Double_t GeneralValues::GetVal(Int_t index)
{
    return fitValues[index];
}

#endif