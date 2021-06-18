/*
    PURPOSE: used to store the single element fitted histograms for a cycle
*/
#ifndef SingleCycleHistoHolder_H
#define SingleCycleHistoHolder_H

#include "TH1D.h"
#include "CycleHistoHolder.h"

using namespace std;

class SingleCycleHistoHolder{
private:
    Int_t numCycles, numRuns, numBins, timeRunEnd, numElements;
    CycleHistoHolder** histoArr;
    string histoName;
    string* elementNames;
public:
    SingleCycleHistoHolder(Int_t numCycles, Int_t numElements, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd, Int_t individualFitChoice, string* elementNames);
    ~SingleCycleHistoHolder();
    TH1D* GetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex);
    void SetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex, TH1D* histo);
};

SingleCycleHistoHolder::SingleCycleHistoHolder(Int_t numCycles, Int_t numElements, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd, Int_t individualFitChoice, string* elementNames)
{
    this->numCycles = numCycles;
    this->numElements = numElements;
    this->numRuns = numRuns;
    this->histoName = histoName;
    this->numBins = numBins;
    this->timeRunEnd = timeRunEnd;
    this->elementNames = elementNames;
    string histoParameter = histoName;
    
    histoArr = new CycleHistoHolder* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        histoName = elementNames[i] + " " + histoName;
        histoArr[i] = new CycleHistoHolder(numCycles, numRuns, histoName, numBins, timeRunEnd, individualFitChoice);
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