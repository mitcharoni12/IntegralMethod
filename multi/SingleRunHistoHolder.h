/*
    PURPOSE: used to store the single element fitted histograms for a run
    looking at the execution of the program this class is acctually never used
*/
#ifndef SINGLERUNHISTOHOLDER_H
#define SINGLERUNHISTOHOLDER_H

#include "TH1D.h"
#include "RunHistoHolder.h"

using namespace std;

class SingleRunHistoHolder{
private:
    Int_t cycleNum, numElements, numRuns, numBins, timeRunEnd;
    RunHistoHolder** singleHistos;
    string histoName;
    string* elementNames;
public:
    SingleRunHistoHolder(Int_t cycleNum, Int_t numElements, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd, string* elementNames);
    ~SingleRunHistoHolder();
    TH1D* GetAHisto(Int_t histoIndex, Int_t elementIndex);
    void SetAHisto(Int_t histoIndex, Int_t elementIndex, TH1D* histo);
};

SingleRunHistoHolder::SingleRunHistoHolder(Int_t cycleNum, Int_t numElements, Int_t numRuns, string histoName, Int_t numBins, Int_t timeRunEnd, string* elementNames)
{
    this->cycleNum = cycleNum;
    this->numElements = numElements;
    this->numRuns = numRuns;
    this->histoName = histoName;
    this->numBins = numBins;
    this->timeRunEnd = timeRunEnd;
    this->elementNames = elementNames;
    string histoParameter = histoName;
    
    singleHistos = new RunHistoHolder* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        histoName = elementNames[i] + " " + histoName;
        cout << "HERE " << histoName << endl << endl << endl << endl;
        singleHistos[i] = new RunHistoHolder(cycleNum, numRuns, histoName, numBins, timeRunEnd);
        histoName = histoParameter;
    }
}

SingleRunHistoHolder::~SingleRunHistoHolder()
{
    delete [] singleHistos;
}

TH1D* SingleRunHistoHolder::GetAHisto(Int_t histoIndex, Int_t elementIndex)
{
    return singleHistos[elementIndex]->GetAHisto(histoIndex);
}

void SingleRunHistoHolder::SetAHisto(Int_t histoIndex, Int_t elementIndex, TH1D* histo)
{
    singleHistos[elementIndex]->SetAHisto(histoIndex, histo);
}

#endif