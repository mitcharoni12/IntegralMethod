#ifndef RUNINTEGRALHISTOHOLDER_H
#define RUNINTEGRALHISTOHOLDER_H

#include "TH1D.h"

using namespace std;

/// Used to store all the histograms for the multiple integral fits in a single cycle for a single element. These histograms are the ones filled with simulated data then fitted.
/// EX: If we are doing 20 cycles and 20 runs per cycles of the decay chain 144Cs->144Ba->144La. This would hold the 20 histograms of a single cycle for the element La.
class RunIntegralHistoHolder{
private:
    TH1D** histoArr;
public:
    RunIntegralHistoHolder(Int_t cycleNum, Int_t numRuns, string histoName, Int_t numBins, Double_t* binEdges);
    ~RunIntegralHistoHolder();
    TH1D* GetAHisto(Int_t histoIndex);
    void SetAHisto(Int_t runIndex, TH1D* histo);
};

/// Constructor for class. Dynaically allocates histograms and names them properly.
RunIntegralHistoHolder::RunIntegralHistoHolder(Int_t cycleNum, Int_t numRuns, string histoName, Int_t numBins, Double_t* binEdges)
{
    string tempHistoName;

    histoArr = new TH1D* [numRuns];
    for(int i = 0; i < numRuns; i++)
    {
        tempHistoName = histoName + " " + "Run: " + (to_string(i+1)).c_str() + " Cycle: " + (to_string(cycleNum+1)).c_str();
        //tempHistoName = "Cumulative Sum Histogram";
        histoArr[i] = new TH1D(tempHistoName.c_str(), tempHistoName.c_str(), numBins, binEdges);
        histoArr[i]->GetXaxis()->SetTitle("Time(S)");
        histoArr[i]->GetYaxis()->SetTitle("Counts");
    }
}

RunIntegralHistoHolder::~RunIntegralHistoHolder()
{
    delete [] histoArr;
}

/// Gets a histogram at a specific run index.
TH1D* RunIntegralHistoHolder::GetAHisto(Int_t runIndex)
{
    TH1D* tempHisto = histoArr[runIndex];
    return tempHisto;
}

/// Sets a histogram at a specific run index.
void RunIntegralHistoHolder::SetAHisto(Int_t runIndex, TH1D* histo)
{
    histoArr[runIndex] = histo;
}

#endif