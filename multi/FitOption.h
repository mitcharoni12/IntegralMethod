#ifndef FITOPTION_H
#define FITOPTION_H

#include "TMath.h"

using namespace std;

/// Used to hold all the fit options for the program. Holds these options reguarless of program execution type chosen.
class FitOption{
private:
    Int_t numRuns = 1, numCycles = 1, numElements, numBins, binRebinInc = 0, startBinLeaveOut = 0, endBinLeaveOut = 0, programExecutionType = 1
         ,eventDecrement = 0, timeShiftType = 1, rebinDifference = 0, binTimeFitInc = 0;
    Double_t timeRunEnd, events, binWidth;
    Double_t* timeFitEndArr, *timeFitStartArr, *binWidthArr;
    string* elementNames;
    Int_t* binNumArr;
    bool multiSource = false, rebin = false, runMeanDifference = false;
public:
    //getters and setters
    void SetNumRuns(Int_t numRuns){this->numRuns = numRuns;}
    void SetNumCycles(Int_t numCycles){this->numCycles = numCycles;}
    void SetNumElements(Int_t numElements){this->numElements = numElements;}
    void SetNumBins(Int_t numBins){this->numBins = numBins;}
    void SetRebinBinInc(Int_t binRebinInc){this->binRebinInc = binRebinInc;}
    void SetLeaveOutStartBinNumber(Int_t startBinLeaveOut){this->startBinLeaveOut = startBinLeaveOut;}
    void SetLeaveOutEndBinNumber(Int_t endBinLeaveOut){this->endBinLeaveOut = endBinLeaveOut;}
    void SetProgramExecutionType(Int_t programExecutionType){this->programExecutionType = programExecutionType;}
    void SetEventDecrement(Int_t eventDecrement){this->eventDecrement = eventDecrement;}
    void SetTimeShiftType(Int_t timeShiftType){this->timeShiftType = timeShiftType;}
    void SetRebinDifference(Int_t rebinDifference){this->rebinDifference = rebinDifference;}
    void SetTimeRunEnd(Double_t timeRunEnd){this->timeRunEnd = timeRunEnd;}
    void SetTimeFitBinInc(Double_t binTimeFitInc){this->binTimeFitInc = binTimeFitInc;}
    void SetNumEvents(Double_t events){this->events = events;}
    void SetBinWidth(Double_t binWidth){this->binWidth = binWidth;}
    void SetMultiSourceChoice(bool multiSource){this->multiSource = multiSource;}
    void SetRebinChoice(bool rebin){this->rebin = rebin;}
    void SetRunMeanDifference(bool runMeanDifference){this->runMeanDifference = runMeanDifference;cout << "SETTING TRUE" << endl;}
    void SetElementNames(string* elementNames){this->elementNames = elementNames;}
    Int_t GetNumRuns(){return numRuns;}
    Int_t GetNumCycles(){return numCycles;}
    Int_t GetNumElements(){return numElements;}
    Int_t GetNumBins(){return numBins;}
    Int_t GetRebinBinInc(){return binRebinInc;}
    Int_t GetStartBinLeaveOut(){return startBinLeaveOut;}
    Int_t GetEndBinLeaveOut(){return endBinLeaveOut;}
    Int_t GetProgramExecutionType(){return programExecutionType;}
    Int_t GetEventDecrement(){return eventDecrement;}
    Int_t GetTimeShiftType(){return timeShiftType;}
    Int_t GetRebinDifference(){return rebinDifference;}
    Int_t* GetBinNumArr(){return binNumArr;}
    Double_t GetTimeRunEnd(){return timeRunEnd;}
    Double_t GetTimeFitBinInc(){return binTimeFitInc;}
    Double_t GetNumEvents(){return events;}
    Double_t GetBinWidth(){return binWidth;}
    Double_t* GetTimeFitEndArr(){return timeFitEndArr;}
    Double_t* GetTimeFitStartArr(){return timeFitStartArr;}
    Double_t* GetBinWidthArr(){return binWidthArr;}
    bool GetMultiSourceChoice(){return multiSource;}
    bool GetRebinChoice(){return rebin;}
    bool GetRunMeanDifference(){return runMeanDifference;}
    string* GetElementNames(){return elementNames;}
    //functions
    void CreateRequiredDataSets();
    ~FitOption();
};

/// Creates any required data set such as bin number array, time fit end array, time fit start array, and bin width array.
/// These data sets are needed to be generated so that the things like the histograms can be generated.
void FitOption::CreateRequiredDataSets()
{
    binNumArr = new Int_t [numCycles];
    timeFitEndArr = new Double_t [numCycles];
    binWidthArr = new Double_t [numCycles];
    timeFitStartArr = new Double_t [numCycles];

    //need to calculate bin width initally
    if(rebin)
    {
        Double_t rebinSize = (Double_t) numBins;
        Double_t tempWidth;

        for(int i = 0; i < numCycles; i++)
        {
            timeFitEndArr[i] = timeRunEnd;
            tempWidth = timeRunEnd / rebinSize;
            binWidthArr[i] = tempWidth;
            binNumArr[i] = rebinSize;
            timeFitStartArr[i] = 0.0f;
            rebinSize = rebinSize + rebinDifference;
            tempWidth = timeRunEnd / rebinSize;
        }
    }
    //need to calculate number of bins initally, change time fit based on number bins
    if(!rebin)
    {
        Double_t fitEnd = timeRunEnd;
        Double_t fitStart = 0.0f;
        Double_t initBinNum = fitEnd / binWidth;
        Double_t addedBins = 0.0f;
        Double_t totalBins = initBinNum;
        for(int i = 0; i < numCycles; i++)
        {
            timeFitEndArr[i] = fitEnd;
            binNumArr[i] = totalBins;
            binWidthArr[i] = binWidth;
            timeFitStartArr[i] = fitStart;

            cout << "FIT END: " << fitEnd << " NUM BINS: " << totalBins << " BIN WIDTH: " << binWidth << " START FIT: " << fitStart << endl;
            addedBins = addedBins + binTimeFitInc;
            //move end
            if(timeShiftType == 1){
                totalBins = initBinNum + addedBins;
                fitEnd = totalBins * binWidth;
            //move start
            }else if(timeShiftType == 2)
            {
                fitStart = addedBins * binWidth;
            }
        }
    }
}

FitOption::~FitOption()
{
    delete [] binNumArr;
    delete [] binWidthArr;
    delete [] timeFitEndArr;
}

#endif