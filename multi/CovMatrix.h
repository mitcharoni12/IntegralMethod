#ifndef COVMATRIX_H
#define COVMATRIX_H

#include "TMath.h"
#include "TMatrixD.h"

using namespace std;

class CovMatrix{
private:
    Double_t* batemanHistoCov;
    Double_t* integralHistoCov;
    TMatrixD* covMatrix;
    Int_t numBins;
public:
    CovMatrix(Int_t numBins);
    ~CovMatrix();

};

CovMatrix::CovMatrix(Int_t numBins, Double_t* batemanHistoCov, Double_t* integralHistoCov)
{
    this->batemanHistoCov = batemanHistoCov;
    this->integralHistoCov = integralHistoCov;
    this->numBins = numBins;
    
}

#endif