#include <iostream>
#include "TF1.h"

using namespace std;

typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);

Double_t testFunction(Double_t *x, Double_t* par);

void compile()
{
    decayFunction* regular = new decayFunction [1];
    regular[0] = testFunction;
    TF1* testyBoi = new TF1("boi", regular[0], 0., 10.);
    testyBoi->Draw();
    delete [] regular;
}

Double_t testFunction(Double_t *x, Double_t* par)
{
    Double_t time = x[0];
    Double_t f = time;
    return f;
}
