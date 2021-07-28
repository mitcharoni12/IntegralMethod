#ifndef SINGLECYCLEHISTOHOLDER_H
#define SINGLECYCLEHISTOHOLDER_H

#include "TH1D.h"
#include "CycleHistoHolder.h"

using namespace std;

/// Holds all the histograms for all the single histograms of the program.
/// EX: If we are doing 20 cycles of 20 runs for the decay chain 144Cs->144Ba->144La we have single histograms for each element in the decay chain and a total integral and bateman histogram.
/// The total integral and Bateman histograms for the 20 cycles of 20 runs are held by the CylceHistoHolder class but the single histograms for Cs, Ba, and La are all held by this class.
class SingleCycleHistoHolder{
private:
    Int_t numElements;
    CycleHistoHolder** histoArr; ///< array of size numElements of CycleHistoHolder objects.
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

/// Gets a histogram for a specific run, cycle, and element index.
TH1D* SingleCycleHistoHolder::GetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex)
{
    return histoArr[elementIndex]->GetAHisto(cycleIndex, runIndex);
}

/// Sets a histogram for a specific run, cycle, and element index.
void SingleCycleHistoHolder::SetAHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex, TH1D* histo)
{
    histoArr[elementIndex]->SetAHisto(cycleIndex, runIndex, histo);
}

#endif