#ifndef FITPARAMETERSTORE_H
#define FITPARAMETERSTORE_H

class FitParameterStore{
    public:
        Double_t* getRegularHalfLife(){return regularHalfLife;}
        Double_t* getRegularN0(){return regularN0;}
        Double_t* getIntegralHalfLife(){return integralHalfLife;}
        Double_t* getIntegralN0(){return integralN0;}
        Double_t* getRegularHalfLifeError(){return regularHalfLifeError;}
        Double_t* getRegularN0Error(){return regularN0Error;}
        Double_t* getIntegralHalfLifeError(){return integralHalfLifeError;}
        Double_t* getIntegralN0Error(){return integralN0Error;}
        FitParameterStore(Int_t numElements_p);
        ~FitParameterStore();
    private:
        Int_t numElements;
        Double_t *regularHalfLife;
        Double_t *regularN0;
        Double_t *integralHalfLife;
        Double_t *integralN0;
        Double_t *regularHalfLifeError;
        Double_t *regularN0Error;
        Double_t *integralHalfLifeError;
        Double_t *integralN0Error;
};

FitParameterStore::FitParameterStore(Int_t numElements_p)
{
    this->numElements = numElements_p;
    regularHalfLife = new Double_t [numElements];
    regularN0 = new Double_t [numElements];
    integralHalfLife = new Double_t [numElements];
    integralN0 = new Double_t [numElements];
    regularHalfLifeError = new Double_t [numElements];
    regularN0Error = new Double_t [numElements];
    integralHalfLifeError = new Double_t [numElements];
    integralN0Error = new Double_t [numElements];
}

FitParameterStore::~FitParameterStore()
{
    delete regularHalfLife;
    delete regularN0;
    delete integralHalfLife;
    delete integralN0;
    delete regularHalfLifeError;
    delete regularN0Error;
    delete integralHalfLifeError;
    delete integralN0Error;
}

#endif