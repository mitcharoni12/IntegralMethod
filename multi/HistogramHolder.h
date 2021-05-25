/*
    PURPOSE: used to store the fitted histograms for each fit in a run or cycle
*/
#ifndef HISTOGRAMHOLDER_H
#define HISTOGRAMHOLDER_h

#include "TH1D.h"

using namespace std;

class HistogramHolder{
private:
    TH1D* histoArr;
    Int_t numFits, numBins, timeRunEnd;
    string histoName;
public:
    HistogramHolder(Int_t numFits, string histoName, Int_t numBins, Int_t timeRunEnd);
    ~HistogramHolder();
};

HistogramHolder::HistogramHolder(Int_t numFits, string histoName, Int_t numBins, Int_t timeRunEnd)
{
    this->numFits = numFits;
    this->numBins = numBins;
    this->numTimeEnd = timeRunEnd;
    this->histoName = histoName;
    histoArr = new TH1D* [numFits];
    for(int i = 0; i < numFits; i++)
    {
        histoName = histoName + " " + (i+1);
        histoArr[i] = new TH1D(histoName.c_str(), histoName.c_str(), numBins, 0., timeRunEnd);
    }
}

HistogramHolder::~HistogramHolder()
{
    delete [] histoArr;
}

#endif