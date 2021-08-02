#ifndef FITOPTION_H
#define FITOPTION_H

#include <math.h>
#include "TMath.h"

using namespace std;

/// Used to hold all the fit options for the program. Holds these options reguarless of program execution type chosen.
class FitOption{
private:
    Int_t numRuns = 1, numCycles = 1, numElements, numBins, binRebinInc = 0, leaveOutStartBinsSim = 0, leaveOutEndBinsSim = 0, programExecutionType = 1,
          leaveOutStartBinsInput = 0, leaveOutEndBinsInput = 0,eventDecrement = 0, timeShiftType = 1, rebinDifference = 0, binTimeFitInc = 0, inputHistoExecutionType = 1
          ,inputHistoBinNum;
    Double_t timeRunEndSimulated, timeRunStartInput, timeRunEndInput, events, binWidth, inputHistoTimeEnd, inputTimeInc;
    Double_t* timeFitEndArr, *timeFitStartArr, *binWidthArr, *timeLengthArr, *fitStartBinArr, *fitEndBinArr;
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
    void SetLeaveOutStartBinsSim(Int_t leaveOutStartBinsSim){this->leaveOutStartBinsSim = leaveOutStartBinsSim;}
    void SetLeaveOutEndBinsSim(Int_t leaveOutEndBinsSim){this->leaveOutEndBinsSim = leaveOutEndBinsSim;}
    void SetLeaveOutStartBinsInput(Int_t leaveOutStartBinsInput){this->leaveOutStartBinsInput = leaveOutStartBinsInput;}
    void SetLeaveOutEndBinsInput(Int_t leaveOutEndBinsInput){this->leaveOutEndBinsInput = leaveOutEndBinsInput;}
    void SetProgramExecutionType(Int_t programExecutionType){this->programExecutionType = programExecutionType;}
    void SetEventDecrement(Int_t eventDecrement){this->eventDecrement = eventDecrement;}
    void SetTimeShiftType(Int_t timeShiftType){this->timeShiftType = timeShiftType;}
    void SetRebinDifference(Int_t rebinDifference){this->rebinDifference = rebinDifference;}
    void SetInputHistoExecutionType(Int_t inputHistoExecutionType){this->inputHistoExecutionType = inputHistoExecutionType;}
    void SetInputHistoBinNum(Int_t inputHistoBinNum){this->inputHistoBinNum = inputHistoBinNum;}
    void SetTimeRunEndSimulated(Double_t timeRunEndSimulated){this->timeRunEndSimulated = timeRunEndSimulated;}
    void SetTimeRunEndInput(Double_t timeRunEndInput){this->timeRunEndInput = timeRunEndInput;}
    void SetTimeRunStartInput(Double_t timeRunStartInput){this->timeRunStartInput = timeRunStartInput;}
    void SetTimeFitBinInc(Double_t binTimeFitInc){this->binTimeFitInc = binTimeFitInc;}
    void SetNumEvents(Double_t events){this->events = events;}
    void SetBinWidth(Double_t binWidth){this->binWidth = binWidth;}
    void SetInputHistoTimeEnd(Double_t inputHistoTimeEnd){this->inputHistoTimeEnd = inputHistoTimeEnd;}
    void SetInputTimeInc(Double_t inputTimeInc){this->inputTimeInc = inputTimeInc;}
    void SetMultiSourceChoice(bool multiSource){this->multiSource = multiSource;}
    void SetRebinChoice(bool rebin){this->rebin = rebin;}
    void SetRunMeanDifference(bool runMeanDifference){this->runMeanDifference = runMeanDifference;}
    void SetElementNames(string* elementNames){this->elementNames = elementNames;}
    Int_t GetNumRuns(){return numRuns;}
    Int_t GetNumCycles(){return numCycles;}
    Int_t GetNumElements(){return numElements;}
    Int_t GetNumBins(){return numBins;}
    Int_t GetRebinBinInc(){return binRebinInc;}
    Int_t GetLeaveOutStartBinsSim(){return leaveOutStartBinsSim;}
    Int_t GetLeaveOutEndBinsSim(){return leaveOutEndBinsSim;}
    Int_t GetLeaveOutStartBinsInput(){return leaveOutStartBinsInput;}
    Int_t GetLeaveOutEndBinsInput(){return leaveOutEndBinsInput;}
    Int_t GetProgramExecutionType(){return programExecutionType;}
    Int_t GetEventDecrement(){return eventDecrement;}
    Int_t GetTimeShiftType(){return timeShiftType;}
    Int_t GetRebinDifference(){return rebinDifference;}
    Int_t GetInputHistoExecutionType(){return inputHistoExecutionType;}
    Int_t GetInputHistoBinNum(){return inputHistoBinNum;}
    Int_t* GetBinNumArr(){return binNumArr;}
    Double_t GetTimeRunEndSimulated(){return timeRunEndSimulated;}
    Double_t GetTimeRunEndInput(){return timeRunEndInput;}
    Double_t GetTimeRunStartInput(){return timeRunStartInput;}
    Double_t GetTimeFitBinInc(){return binTimeFitInc;}
    Double_t GetNumEvents(){return events;}
    Double_t GetBinWidth(){return binWidth;}
    Double_t GetInputHistoTimeEnd(){return inputHistoTimeEnd;}
    Double_t GetInputTimeInc(){return inputTimeInc;}
    Double_t* GetTimeFitEndArr(){return timeFitEndArr;}
    Double_t* GetTimeFitStartArr(){return timeFitStartArr;}
    Double_t* GetBinWidthArr(){return binWidthArr;}
    Double_t* GetTimeLengthArr(){return timeLengthArr;}
    Double_t* GetFitStartBinArr(){return fitStartBinArr;}
    Double_t* GetFitEndBinArr(){return fitEndBinArr;}
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
    timeLengthArr = new Double_t [numCycles];
    fitStartBinArr = new Double_t [numCycles];
    fitEndBinArr = new Double_t [numCycles];

    //need to calculate bin width initally
    if(rebin)
    {
        Double_t rebinSize = (Double_t) numBins;
        Double_t tempWidth;

        for(int i = 0; i < numCycles; i++)
        {
            timeFitEndArr[i] = timeRunEndSimulated;
            tempWidth = timeRunEndSimulated / rebinSize;
            binWidthArr[i] = tempWidth;
            binNumArr[i] = rebinSize;
            timeFitStartArr[i] = 0.0f;
            rebinSize = rebinSize + rebinDifference;
            tempWidth = timeRunEndSimulated / rebinSize;
        }
    }
    //need to calculate number of bins initally, change time fit based on number bins
    if(!rebin && inputHistoExecutionType == 1)
    {
        Double_t fitEnd = timeRunEndSimulated;
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
    if(inputHistoExecutionType == 3 || inputHistoExecutionType == 2)
    {
        Double_t binWidth = inputHistoTimeEnd / inputHistoBinNum;
        Double_t fitStart = timeRunStartInput;
        Double_t fitEnd = timeRunEndInput;
        Double_t fitStartBinNum;
        Double_t fitEndBinNum;
        for(int i = 0; i < numCycles; i++)
        {
            timeFitEndArr[i] = fitEnd;
            timeFitStartArr[i] = fitStart;
            timeLengthArr[i] = fitEnd - fitStart;
            binWidthArr[i] = binWidth;
            //I determined this was the most efficient way to get number of bins when cutting out the input histogram.
            //We need this becuase we cut the input histogram and put it in a new histogram to make
            fitStartBinNum = floor(fitStart / binWidth);
            fitEndBinNum = floor(fitEnd / binWidth);
            fitStartBinArr[i] = fitStartBinNum;
            fitEndBinArr[i] = fitEndBinNum;
            binNumArr[i] = fitEndBinNum - fitStartBinNum + 1;

            cout << "FIT END: " << fitEnd << " NUM BINS: " << binNumArr[i] << " BIN WIDTH: " << binWidth << " START FIT: " << fitStart << endl;
            //move end
            if(timeShiftType == 1)
            {
                fitEnd = fitEnd + inputTimeInc;
            //move start
            }else if(timeShiftType == 2)
            {
                fitStart = fitStart + inputTimeInc;
            }
        }
    }
}

FitOption::~FitOption()
{
    delete [] binNumArr;
    delete [] binWidthArr;
    delete [] timeFitEndArr;
    delete [] timeFitStartArr;
    delete [] timeLengthArr;
    delete [] fitStartBinArr;
    delete [] fitEndBinArr;
}

#endif