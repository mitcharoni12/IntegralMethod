#ifndef GENERALVALUES_H
#define GENERALVALUES_H

#include "TMath.h"

/// Like a glorified array, purpose is to store values extracted from multiple runs.
/// Ex: if we do 20x runs, we want a way to store each half life in an array life format. To graph the fitted values for the 20x runs we have to have it in this array life format.
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

/// sets an array value
void GeneralValues::SetVal(Int_t index, Double_t val)
{
    fitValues[index] = val;
}

/// gets an array value
Double_t GeneralValues::GetVal(Int_t index)
{
    return fitValues[index];
}

#endif