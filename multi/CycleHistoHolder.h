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
    Int_t numCycles, numRuns, timeRunEnd;
    Int_t* binArr
    string histoName;
    RunHistoHolder** histoArr;
public:
    CycleHistoHolder(Int_t numCycles, Int_t numRuns, string histoName, Int_t* binArr, Int_t timeRunEnd, Int_t individualFitChoice);
    ~CycleHistoHolder();
    TH1D* GetAHisto(Int_t cycleIndex, Int_t runIndex);
    void SetAHisto(Int_t cycleIndex, Int_t runIndex, TH1D* histo);
};

CycleHistoHolder::CycleHistoHolder(Int_t numCycles, Int_t numRuns, string histoName, Int_t* binArr, Int_t timeRunEnd, Int_t individualFitChoice)
{
    this->numCycles = numCycles;
    this->numRuns = numRuns;
    this->binArr = binArr;
    this->timeRunEnd = timeRunEnd;
    this->histoName = histoName;

    histoArr = new RunHistoHolder* [numCycles];
    for(int i = 0; i < numCycles; i++)
    {
        histoArr[i] = new RunHistoHolder(i, numRuns, histoName, binArr[i], timeRunEnd, individualFitChoice);
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