/*
    PURPOSE: used to group together some of the fit options used in the program
*/
#ifndef FITOPTION_H
#define FITOPTION_H

#include "TMath.h"

class FitOption{
private:
    Int_t numRuns = 1, numCycles = 1, numElements, numBins, binRebinInc = 0, startBinLeaveOut = 0, endBinLeaveOut = 0, programExecutionType = 1
         ,eventDecrement = 0, timeShiftType, rebinDifference = 0, binTimeFitInc = 0;
    Double_t timeRunEnd, events, binWidth;
    Double_t* timeFitEndArr;
    Int_t* rebinBinNumbers;
    bool multiSource = false, rebin = false, runMeanDifference = false;
public:
    void SetNumRuns(Int_t numRuns){this->numRuns = numRuns;}
    void SetNumCycles(Int_t numCycles){this->numCycles = numCycles;}
    void SetNumElements(Int_t numElements){this->numElements = numElements;}
    void SetNumEvents(Int_t numEvents){this->numEvents = numEvents;}
    void SetRebinBinInc(Int_t binRebinInc){this->binRebinInc = binRebinInc;}
    void SetStartBinLeaveOut(Int_t startBinLeaveOut){this->startBinLeaveOut = startBinLeaveOut;}
    void SetEndBinLeaveOut(Int_t endBinLeaveOut){this->endBinLeaveOut = endBinLeaveOut;}
    void SetProgramExecutionType(Int_t SetProgramExecutionType){this->programExecutionType = programExecutionType;}
    void SetEventDecrement(Int_t eventDecrement){this->eventDecrement = eventDecrement;}
    void SetTimeShiftType(Int_t timeShiftType){this->timeShiftType = timeShiftType;}
    void SetRebinDifference(Int_t rebinDifference){this->rebinDifference = rebinDifference;}
    void SetTimeRunEnd(Double_t timeRunEnd){this->timeRunEnd = timeRunEnd;}
    void SetTimeFitBinInc(Double_t binTimeFitInc){this->binTimeFitInc = binTimeFitInc;}
    void SetEvents(Double_t events){this->events = events;}
    void SetBinWidth(Double_t binWidth){this->binWidth = binWidth;}
    void SetMultiSourceChoice(bool multiSource){this->multiSource = multiSource;}
    void SetRebinChoice(bool rebin){this->rebin = rebin;}
    void SetRunMeanDifference(bool runMeanDifference){this->runMeanDifference = runMeanDifference;}
    Int_t GetNumRuns(){return numRuns;}
    Int_t GetNumEvents(){return numElements;}
    Int_t GetNumBins(){return numBins;}
    Int_t GetNumEvents(){return numEvents;}
    Int_t GetBinInc(){return binInc;}
    Int_t GetStartBinLeaveOut(){return startBinLeaveOut;}
    Int_t GetEndBinLeaveOut(){return endBinLeaveOut;}
    Int_t GetProgramExecutionType(){return programExecutionType;}
    Int_t GetEventDecrement(){return eventDecrement;}
    Int_t GetTimeShiftType(){return timeShiftType;}
    Int_t GetRebinDifference(){return rebinDifference;}
    Double_t GetTimeRunEnd(){return timeRunEnd;}
    Double_t GetTimeInc(){return timeInc;}
    Double_t GetEvents(){return events;}
    Double_t GetBinWidth(){return binWidth;}
    bool GetMultiSourceChoice(){return multiSource;}
    bool GetRebinChoice(){return rebin;}
    bool GetRunMeanDifference(){return runMeanDifference;}
};

FitOptions::CreateRequiredDataSets()
{
    rebinBinNumbers = new Int_t [numCycles];

    Int_t rebinSize = numBins;
    for(int i = 0; i < numCycles; i++)
    {
        rebinBinNumbers[i] = rebinSize;
        rebinSize = rebinSize + rebinDifference;
    }



}

FitOptions::~FitOptions()
{
    delete [] rebinBinNumbers;
}

#endif