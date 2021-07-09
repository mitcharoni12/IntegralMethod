/*
    PURPOSE: used to group together some of the fit options used in the program
*/
#ifndef FITOPTIONS_H
#define FITOPTIONS_H

#include "TMath.h"

class FitOptions{
private:
    Int_t numRuns, numCycles, numElements, numBins, binInc, startBinLeaveOut, endBinLeaveOut;
    Double_t timeRunEnd, timeInc, events;
    bool multiSource; rebin;
public:
    void SetNumRuns(Int_t numRuns){this->numRuns = numRuns;}
    void SetNumCycles(Int_t numCycles){this->numCycles = numCycles;}
    void SetNumElements(Int_t numElements){this->numElements = numElements;}
    void SetBinInc(Int_t binInc){this->binInc = binInc;}
    void SetStartBinLeaveOut(Int_t startBinLeaveOut){this->startBinLeaveOut = startBinLeaveOut;}
    void SetEndBinLeaveOut(Int_t endBinLeaveOut){this->endBinLeaveOut = endBinLeaveOut;}
    void SetTimeRunEnd(Double_t timeRunEnd){this->timeRunEnd = timeRunEnd;}
    void SetTimeInc(Double_t timeInc){this->timeInc = timeInc;}
    void SetEvents(Double_t events){this->events = events;}
    void SetMultiSourceChoice(bool multiSource){this->multiSource = multiSource;}
    void SetRebinChoice(bool rebin){this->rebin = rebin;}
};



#endif