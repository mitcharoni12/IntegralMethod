#ifndef CYCLEINTEGRALHISTOHOLDER_H
#define CYCLEINTEGRALHISTOHOLDER_H

using namespace std;

#include "TH1D.h"
#include "RunIntegralHistoHolder.h"

/// Used to store the histograms for the multiple integral cycles of a single element. These histograms contain the simulated data and the ones that are fit.
/// EX: If we are doing 20 cycles of 20 runs for the decay chain 144Cs->144Ba->144La this class would hold all the histograms for the 20 cycles of La, each cycle held by a RunIntegralHistoHolder class.
class CycleHistoHolder{
private:
    Int_t numCycles;
    RunIntegralHistoHolder** histoArr; ///< Holds all the histograms, array of RunIntegralHistoHolder objects.
public:
    CycleHistoHolder(Int_t numCycles, Int_t numRuns, string histoName, Int_t* binArr, Double_t** binEdgesArr);
    ~CycleHistoHolder();
    TH1D* GetAHisto(Int_t cycleIndex, Int_t runIndex);
    void SetAHisto(Int_t cycleIndex, Int_t runIndex, TH1D* histo);
};

CycleHistoHolder::CycleHistoHolder(Int_t numCycles, Int_t numRuns, string histoName, Int_t* binArr, Double_t** binEdgesArr)
{
    this->numCycles = numCycles;

    histoArr = new RunIntegralHistoHolder* [numCycles];
    for(int i = 0; i < numCycles; i++)
    {
        histoArr[i] = new RunIntegralHistoHolder(i, numRuns, histoName, binArr[i], binEdgesArr[i]);
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

/// Gets a histogram at a specific run and cycle index.
TH1D* CycleHistoHolder::GetAHisto(Int_t cycleIndex, Int_t runIndex)
{
    return histoArr[cycleIndex]->GetAHisto(runIndex);
}

/// Sets a histogram at a specific run and cycle index.
void CycleHistoHolder::SetAHisto(Int_t cycleIndex, Int_t runIndex, TH1D* histo)
{
    histoArr[cycleIndex]->SetAHisto(runIndex, histo);
}

#endif