#ifndef FITVALUES_H
#define FITVALUES_H

#include "TMath.h"

/// Used for storing fitted values for a single element componet in a decay chain.
/// EX: Decay chain 144Cs->144Ba->144La, each element componet in that chain have parameters values for N0, N0Error, half life, and half life error.
/// Lets look at the values extracted in the single La histogram, the bateman equation for La also has parameters for N0, and half life of Cs and Ba.
/// So this class would for example store the N0, N0Error, half life, and half life error of Ba in the single La histogram.
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