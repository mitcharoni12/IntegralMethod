/*
    PURPOSE: used to store the fitted histograms for each fit in a run
*/
#ifndef RUNHISTOHOLDER_H
#define RUNHISTOHOLDER_H

#include "TH1D.h"

using namespace std;

class RunHistoHolder{
private:
    TH1D** histoArr;
    Int_t cycleNum, numRuns, numBins, timeRunEnd, individualFitChoice;
    string histoName;
public:
    RunHistoHolder(Int_t cycleNum, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd, Int_t individualFitChoice);
    ~RunHistoHolder();
    TH1D* GetAHisto(Int_t histoIndex);
    void SetAHisto(Int_t runIndex, TH1D* histo);
};

RunHistoHolder::RunHistoHolder(Int_t cycleNum, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd, Int_t individualFitChoice)
{
    string tempHistoName;
    this->individualFitChoice = individualFitChoice;
    this->cycleNum = cycleNum;
    this->numRuns = numRuns;
    this->numBins = numBins;
    this->timeRunEnd = timeRunEnd;
    this->histoName = histoName;

    histoArr = new TH1D* [numRuns];
    for(int i = 0; i < numRuns; i++)
    {
        tempHistoName = histoName + " " + "Run: " + (to_string(i+1)).c_str() + " Cycle: " + (to_string(cycleNum+1)).c_str();
        histoArr[i] = new TH1D(tempHistoName.c_str(), tempHistoName.c_str(), numBins, 0., timeRunEnd);
    }
}

RunHistoHolder::~RunHistoHolder()
{
    if(individualFitChoice != 1)
    {
        for(int i = 0; i < numRuns; i++)
        {
            delete histoArr [i];
        }
    }
    delete [] histoArr;
}

TH1D* RunHistoHolder::GetAHisto(Int_t runIndex)
{
    TH1D* tempHisto = histoArr[runIndex];
    return tempHisto;
}

void RunHistoHolder::SetAHisto(Int_t runIndex, TH1D* histo)
{
    histoArr[runIndex] = histo;
}

#endif