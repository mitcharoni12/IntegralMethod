#include <iostream>
#include "../../../root/root-6.24.00-install/include/TF1.h"

using namespace std;

Double_t testFunction(Double_t *x, Double_t* par);

void compile()
{
    cout << "BRUH";
    TF1* testyBoi = new TF1("boi", testFunction, 0., 10.);
    testyBoi->Draw();   
}

Double_t testFunction(Double_t *x, Double_t* par)
{
    Double_t time = x[0];
    Double_t f = time;
    return f;
}
