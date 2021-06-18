/*
    PURPOSE: used to store the fitted histograms for each fit in a run
*/
#ifndef CYCLEHISTOHOLDER_H
#define CYCLEHISTOHOLDER_H

using namespace std;

#include "TH1D.h"
#include "RunHistoHolder.h"

class CycleHistoHolder{
private:
    Int_t numCycles, numRuns, numBins, timeRunEnd;
    string histoName;
    RunHistoHolder** histoArr;
public:
    CycleHistoHolder(Int_t numCycles, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd);
    ~CycleHistoHolder();
    TH1D* GetAHisto(Int_t cycleIndex, Int_t runIndex);
    void SetAHisto(Int_t cycleIndex, Int_t runIndex, TH1D* histo);
};

CycleHistoHolder::CycleHistoHolder(Int_t numCycles, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd)
{
    this->numCycles = numCycles;
    this->numRuns = numRuns;
    this->numBins = numBins;
    this->timeRunEnd = timeRunEnd;
    this->histoName = histoName;

    histoArr = new RunHistoHolder* [numCycles];
    for(int i = 0; i < numCycles; i++)
    {
        histoArr[i] = new RunHistoHolder(i, numRuns, histoName, numBins, timeRunEnd);
    }
}

CycleHistoHolder::~CycleHistoHolder()
{
    for(int i = 0; i < numCycles; i++)
    {
        delete histoArr[i];
    }
    delete [] histoArr;
}

TH1D* CycleHistoHolder::GetAHisto(Int_t cycleIndex, Int_t runIndex)
{
    return histoArr[cycleIndex]->GetAHisto(runIndex);
}

void CycleHistoHolder::SetAHisto(Int_t cycleIndex, Int_t runIndex, TH1D* histo)
{
    histoArr[cycleIndex]->SetAHisto(runIndex, histo);
}

#endif