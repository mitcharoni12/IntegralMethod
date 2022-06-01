#ifndef CYCLEHISTOHOLDER_H
#define CYCLEHISTOHOLDER_H

using namespace std;

#include "TH1D.h"
#include "RunHistoHolder.h"

/// Used to store the histograms for the multiple Bateman cycles of a single element. These histograms contain the simulated data and the ones that are fit.
/// EX: If we are doing 20 cycles of 20 runs for the decay chain 144Cs->144Ba->144La this class would hold all the histograms for the 20 cycles of La, each cycle held by a RunHistoHolder class.
class CycleHistoHolder{
private:
    Int_t numCycles;
    RunHistoHolder** histoArr; ///< Holds all the histograms, array of RunHistoHolder objects.
public:
    CycleHistoHolder(Int_t numCycles, Int_t numRuns, string histoName, Int_t* binArr, Double_t* timeEndArr);
    ~CycleHistoHolder();
    TH1D* GetAHisto(Int_t cycleIndex, Int_t runIndex);
    void SetAHisto(Int_t cycleIndex, Int_t runIndex, TH1D* histo);
};

CycleHistoHolder::CycleHistoHolder(Int_t numCycles, Int_t numRuns, string histoName, Int_t* binArr, Double_t* timeEndArr)
{
    this->numCycles = numCycles;

    histoArr = new RunHistoHolder* [numCycles];
    for(int i = 0; i < numCycles; i++)
    {
        histoArr[i] = new RunHistoHolder(i, numRuns, histoName, binArr[i], timeEndArr[i]);
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