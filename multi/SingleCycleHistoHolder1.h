/*
    PURPOSE: used to store the single element fitted histograms for a cycle
*/
#ifndef SingleCycleHistoHolder1_H
#define SingleCycleHistoHolder1_H

#include "TH1D.h"
#include "CycleHistoHolder.h"

using namespace std;

class SingleCycleHistoHolder1{
private:
    Int_t numCycles, numRuns, numBins, timeRunEnd, numElements;
    CycleHistoHolder** histoArr;
    string histoName;
    string* elementNames;
public:
    SingleCycleHistoHolder1(Int_t numCycles, Int_t numElements, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd, string* elementNames);
    ~SingleCycleHistoHolder1();
    TH1D* GetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex);
    void SetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex, TH1D* histo);
};

SingleCycleHistoHolder1::SingleCycleHistoHolder1(Int_t numCycles, Int_t numElements, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd, string* elementNames)
{
    this->numCycles = numCycles;
    this->numElements = numElements;
    this->numRuns = numRuns;
    this->histoName = histoName;
    this->numBins = numBins;
    this->timeRunEnd = timeRunEnd;
    this->elementNames = elementNames;
    
    histoArr = new CycleHistoHolder* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        histoName = elementNames[i] + " " + histoName;
        histoArr[i] = new CycleHistoHolder(numCycles, numRuns, histoName, numBins, timeRunEnd);
    }
}

SingleCycleHistoHolder1::~SingleCycleHistoHolder1()
{
    delete [] histoArr;
}

TH1D* SingleCycleHistoHolder1::GetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex)
{
    return histoArr[elementIndex]->GetAHisto(cycleIndex, runIndex);
}

void SingleCycleHistoHolder1::SetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex, TH1D* histo)
{
    histoArr[elementIndex]->SetAHisto(cycleIndex, runIndex, histo);
}

#endif