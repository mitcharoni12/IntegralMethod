#ifndef RUNHISTOHOLDER_H
#define RUNHISTOHOLDER_H

#include "TH1D.h"

using namespace std;

/// Used to store all the histograms for the multiple Bateman fits in a single cycle for a single element. These histograms are the ones filled with simulated data then fitted.
/// EX: If we are doing 20 cycles and 20 runs per cycles of the decay chain 144Cs->144Ba->144La. This would hold the 20 histograms of a single cycle for the element La.
class RunHistoHolder{
private:
    TH1D** histoArr;
public:
    RunHistoHolder(Int_t cycleNum, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd);
    ~RunHistoHolder();
    TH1D* GetAHisto(Int_t histoIndex);
    void SetAHisto(Int_t runIndex, TH1D* histo);
};

/// Constructor for class. Dynaically allocates histograms and names them properly.
RunHistoHolder::RunHistoHolder(Int_t cycleNum, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd)
{
    string tempHistoName;

    histoArr = new TH1D* [numRuns];
    for(int i = 0; i < numRuns; i++)
    {
        //tempHistoName = histoName + " " + "Run: " + (to_string(i+1)).c_str() + " Cycle: " + (to_string(cycleNum+1)).c_str();
        tempHistoName = "Decay Histogram";
        histoArr[i] = new TH1D(tempHistoName.c_str(), tempHistoName.c_str(), numBins, 0.0, timeRunEnd);
        histoArr[i]->GetXaxis()->SetTitle("Time(S)");
        histoArr[i]->GetYaxis()->SetTitle("Counts");
    }
}

RunHistoHolder::~RunHistoHolder()
{
    delete [] histoArr;
}

/// Gets a histogram at a specific run index.
TH1D* RunHistoHolder::GetAHisto(Int_t runIndex)
{
    TH1D* tempHisto = histoArr[runIndex];
    return tempHisto;
}

/// Sets a histogram at a specific run index.
void RunHistoHolder::SetAHisto(Int_t runIndex, TH1D* histo)
{
    histoArr[runIndex] = histo;
}

#endif