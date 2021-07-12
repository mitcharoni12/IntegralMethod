/*
    PURPOSE: used to store the single element fitted histograms for a cycle
*/
#ifndef SINGLECYCLEHISTOHOLDER_H
#define SINGLECYCLEHISTOHOLDER_H

#include "TH1D.h"
#include "CycleHistoHolder.h"

using namespace std;

class SingleCycleHistoHolder{
private:
    Int_t numElements;
    CycleHistoHolder** histoArr;
    string histoParameter;
public:
    SingleCycleHistoHolder(Int_t numCycles, Int_t numElements, Int_t numRuns, string histoName, Int_t* binArr, Double_t* timeEndArr, string* elementNames);
    ~SingleCycleHistoHolder();
    TH1D* GetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex);
    void SetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex, TH1D* histo);
};

SingleCycleHistoHolder::SingleCycleHistoHolder(Int_t numCycles, Int_t numElements, Int_t numRuns, string histoName, Int_t* binArr, Double_t* timeEndArr, string* elementNames)
{
    this->numElements = numElements;
    histoParameter = histoName;
    histoArr = new CycleHistoHolder* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        histoName = elementNames[i] + " " + histoName;
        histoArr[i] = new CycleHistoHolder(numCycles, numRuns, histoName, binArr, timeEndArr);
        histoName = histoParameter;
    }
}

SingleCycleHistoHolder::~SingleCycleHistoHolder()
{
    for(int i = 0; i < numElements; i++)
    {
        delete histoArr[i];
    }
    delete [] histoArr;
}

TH1D* SingleCycleHistoHolder::GetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex)
{
    return histoArr[elementIndex]->GetAHisto(cycleIndex, runIndex);
}

void SingleCycleHistoHolder::SetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex, TH1D* histo)
{
    histoArr[elementIndex]->SetAHisto(cycleIndex, runIndex, histo);
}

#endif