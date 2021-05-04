#ifndef GENERALVALUES_H
#define GENERALVALUES_H

#include "TMath.h"

class GeneralValues{
private:
    Double_t* fitValues;
public:
    void SetFitValues(Double_t* fitValues){this->fitValues = fitValues;}
    Double_t* GetFitValues(){return fitValues;}
    ~GeneralValues();
};

GeneralValues::~GeneralValues()
{
    delete fitValues;
}

#endif