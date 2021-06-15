/*
    PURPOSE: All Fits have N0, N0Error, halfLife, halfLifeError. Class provides a SINGLE FIT OPERATION FOR ONE ELEMENT.
*/
#ifndef FITVALUES_H
#define FITVALUES_H

#include "TMath.h"

class FitValues{
private:
    Double_t N0;
    Double_t N0Error;
    Double_t halfLife;
    Double_t halfLifeError;
public:
    void SetN0(Double_t N0){this->N0 = N0;}
    void SetN0Error(Double_t N0Error){this->N0Error = N0Error;}
    void SetHalfLife(Double_t halfLife){this->halfLife = halfLife;}
    void SetHalfLifeError(Double_t halfLifeError){this->halfLifeError = halfLifeError;}
    
    Double_t GetN0(){return N0;}
    Double_t GetN0Error(){return N0Error;}
    Double_t GetHalfLife(){return halfLife;}
    Double_t GetHalfLifeError(){return halfLifeError;}
};

#endif