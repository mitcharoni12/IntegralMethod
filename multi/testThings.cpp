#include "FitValues.h"
#include "ChainFitValues.h"
#include "SingleElementFitValues.h"

int testThings()
{
    FitValues* CSFit = new FitValues();
    CSFit->SetN0(1.0);
    CSFit->SetN0Error(2.0);
    CSFit->SetHalfLife(3.0);
    CSFit->SetHalfLifeError(4.0);
    FitValues* BAFit = new FitValues();
    BAFit->SetN0(5.0);
    BAFit->SetN0Error(6.0);
    BAFit->SetHalfLife(7.0);
    BAFit->SetHalfLifeError(8.0);

    ChainFitValues* CSChain = new ChainFitValues(1);
    ChainFitValues* CSBAChain = new ChainFitValues(2);
    FitValues* tempArr [1] = {CSFit};
    FitValues* tempArr2 [2] = {CSFit, BAFit};
    CSChain->SetChainFitValues(tempArr);
    CSBAChain->SetChainFitValues(tempArr2);

    SingleElementFitValues* singles = new SingleElementFitValues(2);
    ChainFitValues* tempArr3 [2]={CSChain, CSBAChain};
    singles->SetElementChainFitValues(tempArr3);


    FitValues* tempFit;
    for(int i = 0; i < singles->GetNumElements(); i++)
    {
        for(int j = 0; j < i+1; j++)
        {
            tempFit = singles->GetChainFitValues(i)->GetSingleElementFitValues(j);
            cout << tempFit->GetN0() << endl;
            cout << tempFit->GetN0Error() << endl;
            cout << tempFit->GetHalfLife() << endl;
            cout << tempFit->GetHalfLifeError() << endl;
        }
    }
    return 0;
}