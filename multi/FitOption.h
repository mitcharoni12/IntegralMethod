#ifndef FITOPTION_H
#define FITOPTION_H

#include <math.h>
#include "TMath.h"

using namespace std;

/// Used to hold all the fit options for the program. Holds these options reguarless of program execution type chosen.
class FitOption{
private:
    Int_t numRuns = 1, numCycles = 1, numElements, numBins, inputHistoBinNum;
    Int_t fitRetries = 1;               ///< Number of retries of the fit before the program moves on
    Int_t rebinBinInc = 0;              ///< Increment of number of bins between cycles for a rebin type program execution.
    Int_t leaveOutStartBinsSim = 0;     ///< Number of bins to leave out at start for data simulation execution type.
    Int_t leaveOutEndBinsSim = 0;       ///< Number of bins to leave out at end for data simulation execution type.
    Int_t programExecutionType = 1;     ///< Execution type for simulation portion of program, 1 = single run, 2 = multiple runs, 3 = multiple cycles.
    Int_t eventChangeType = 1;          ///< If the changing events will change 1 = multapicatively based on factor, 2 = Linearly based on factor.
    Int_t eventChangeFactor = 1;        ///< Factor by which number of events will change.
    Int_t leaveOutStartBinsInput = 0;   ///< Number of bins to leave out at start for data input histogram execution type.
    Int_t leaveOutEndBinsInput = 0;     ///< Number of bins to leave out at end for data input histogram execution type.
    Int_t timeShiftType = 1;            ///< Determines the type of time shift that will happen between cycles of the program, 1 = start: const, end: change, 2 = start: change, end: const.
    Int_t binTimeFitInc = 0;            ///< Bin increment when changing time fit for simulation of the program. Want to keep bin width same to to change time fit add on bin to end.
    Int_t inputHistoExecutionType = 1;  ///< Determines the execution type reguarding the input histogram. 1 = no input histogram, 2 = monte carlo type error evaluation, 3 = input histogram changing time fit.
    Int_t singleElementDataChoice = 1;  ///< Determines if the program will generates and fit single Bateman/integral histograms. 1 = don't generate single histogram, 2 = generate single histograms.
    Double_t timeRunEndSimulated;       ///< Inital time in which data will stop generating/fitting for the simulation execution type of program.(S)
    Double_t timeRunStartInput;         ///< Inital time to start fitting the input histogram(10^-8S).
    Double_t timeRunEndInput;           ///< Inital time to end fitting the input histogram(10^-8S).
    Double_t events;                    ///< Number of events generated for EACH ELEMENT in the decay chain.
    Double_t binWidth;                  ///< Width of the bins for the simulation execution type of the program.
    Double_t inputHistoTimeEnd;         ///< Time in which the input histogram ends(10^-8s).
    Double_t inputTimeInc;              ///< Increment in which time fit will change between cycles of the input histogram.
    Double_t* timeFitEndArr;            ///< Contains the time in which the fit will end between cycles for either the input histogram or the simulated histogram.
    Double_t* timeFitStartArr;          ///< Contains the time in which the fit will start between cycles for either the input histogram or the simulated histogram.
    Double_t* binWidthArr;              ///< Contains the bin width for the histogram in a certain cycle of the program execution for either the input histogram or simulated histograms.
    Double_t* eventNumArr;              ///< Contains the number of events for a histogram in a certain cycle.
    Double_t* timeLengthArr;            ///< Contains the time length of the fit between cycles of the program for the input histogram.
    Double_t* fitStartBinArr;           ///< Contains the bin on which the fit time will start between cycles of the program for the input histogram.
    Double_t* fitEndBinArr;             ///< Contains the bin on which the fit time will end between cycles of the program for the input histogram.
    Double_t** binEdgesArr;             ///< Contains the bin edges for all histograms
    string* elementNames;               ///< Contains the names of the elements in the decay chain.
    Int_t* binNumArr;                   ///< Contains the number of bins for the histogram that is being fitted between cycles of the program for either the input histogram or the simulated histograms.
    bool displayFitAverages = false;    ///< Determines if program will display the average value of the fitted value in the Cycle class.
    bool eventNumChangeChoice = false;  ///< Determines if the number of events will change between cycles.
    bool multiSource = false;           ///< Determines if simulated data is to be just one histogram or multiple histograms. false = one histogram, true =  multiple histograms.
    bool rebinChoice = false;           ///< Determines if the program is to be run with changing the bin number in a set time span(rebinning) or changing the time fit. false = changing time fit, true = rebinning.
    bool runMeanDifference = false;     ///< Determines if the program is going to display the difference between the integral method value and the bateman method value for the fit values. false = display integral and bateman fit values seperatly, true = display difference between values.
    bool timeFitChoice = false;         ///< Determines if the program is going to change the time fit between program cycles
public:
    //getters and setters
    void SetNumRuns(Int_t numRuns){this->numRuns = numRuns;}
    void SetNumCycles(Int_t numCycles){this->numCycles = numCycles;}
    void SetNumElements(Int_t numElements){this->numElements = numElements;}
    void SetNumBins(Int_t numBins){this->numBins = numBins;}
    void SetRebinBinInc(Int_t rebinBinInc){this->rebinBinInc = rebinBinInc;}
    void SetLeaveOutStartBinsSim(Int_t leaveOutStartBinsSim){this->leaveOutStartBinsSim = leaveOutStartBinsSim;}
    void SetLeaveOutEndBinsSim(Int_t leaveOutEndBinsSim){this->leaveOutEndBinsSim = leaveOutEndBinsSim;}
    void SetLeaveOutStartBinsInput(Int_t leaveOutStartBinsInput){this->leaveOutStartBinsInput = leaveOutStartBinsInput;}
    void SetLeaveOutEndBinsInput(Int_t leaveOutEndBinsInput){this->leaveOutEndBinsInput = leaveOutEndBinsInput;}
    void SetProgramExecutionType(Int_t programExecutionType){this->programExecutionType = programExecutionType;}
    void SetEventChangeType(Int_t eventChangeType){this->eventChangeType = eventChangeType;}
    void SetEventChangeFactor(Int_t eventChangeFactor){this->eventChangeFactor = eventChangeFactor;}
    void SetTimeShiftType(Int_t timeShiftType){this->timeShiftType = timeShiftType;}
    void SetInputHistoExecutionType(Int_t inputHistoExecutionType){this->inputHistoExecutionType = inputHistoExecutionType;}
    void SetSingleElementDataChoice(Int_t singleElementDataChoice){this->singleElementDataChoice = singleElementDataChoice;}
    void SetInputHistoBinNum(Int_t inputHistoBinNum){this->inputHistoBinNum = inputHistoBinNum;}
    void SetTimeRunEndSimulated(Double_t timeRunEndSimulated){this->timeRunEndSimulated = timeRunEndSimulated;}
    void SetTimeRunEndInput(Double_t timeRunEndInput){this->timeRunEndInput = timeRunEndInput;}
    void SetTimeRunStartInput(Double_t timeRunStartInput){this->timeRunStartInput = timeRunStartInput;}
    void SetTimeFitBinInc(Double_t binTimeFitInc){this->binTimeFitInc = binTimeFitInc;}
    void SetNumEvents(Double_t events){this->events = events;}
    void SetBinWidth(Double_t binWidth){this->binWidth = binWidth;}
    void SetInputHistoTimeEnd(Double_t inputHistoTimeEnd){this->inputHistoTimeEnd = inputHistoTimeEnd;}
    void SetInputTimeInc(Double_t inputTimeInc){this->inputTimeInc = inputTimeInc;}
    void SetDisplayFitAverages(bool displayFitAverages){this->displayFitAverages = displayFitAverages;}
    void SetEventNumChangeChoice(bool eventNumChangeChoice){this->eventNumChangeChoice = eventNumChangeChoice;}
    void SetMultiSourceChoice(bool multiSource){this->multiSource = multiSource;}
    void SetRebinChoice(bool rebinChoice){this->rebinChoice = rebinChoice;}
    void SetRunMeanDifference(bool runMeanDifference){this->runMeanDifference = runMeanDifference;}
    void SetTimeFitChoice(bool timeFitChoice){this->timeFitChoice = timeFitChoice;}
    void SetElementNames(string* elementNames){this->elementNames = elementNames;}
    Int_t GetNumRuns(){return numRuns;}
    Int_t GetNumCycles(){return numCycles;}
    Int_t GetNumElements(){return numElements;}
    Int_t GetNumBins(){return numBins;}
    Int_t GetRebinBinInc(){return rebinBinInc;}
    Int_t GetLeaveOutStartBinsSim(){return leaveOutStartBinsSim;}
    Int_t GetLeaveOutEndBinsSim(){return leaveOutEndBinsSim;}
    Int_t GetLeaveOutStartBinsInput(){return leaveOutStartBinsInput;}
    Int_t GetLeaveOutEndBinsInput(){return leaveOutEndBinsInput;}
    Int_t GetProgramExecutionType(){return programExecutionType;}
    Int_t GetEventChangeType(){return eventChangeType;}
    Int_t GetEventChangeFactor(){return eventChangeFactor;}
    Int_t GetTimeShiftType(){return timeShiftType;}
    Int_t GetInputHistoExecutionType(){return inputHistoExecutionType;}
    Int_t GetInputHistoBinNum(){return inputHistoBinNum;}
    Int_t GetSingleElementDataChoice(){return singleElementDataChoice;}
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
    Double_t* GetEventNumArr(){return eventNumArr;}
    Double_t* GetTimeLengthArr(){return timeLengthArr;}
    Double_t* GetFitStartBinArr(){return fitStartBinArr;}
    Double_t* GetFitEndBinArr(){return fitEndBinArr;}
    Double_t** GetBinEdges(){return binEdgesArr;}
    bool GetDisplayFitAverages(){return displayFitAverages;}
    bool GetNumEventChangeChoice(){return eventNumChangeChoice;}
    bool GetMultiSourceChoice(){return multiSource;}
    bool GetRebinChoice(){return rebinChoice;}
    bool GetRunMeanDifference(){return runMeanDifference;}
    bool GetTimeFitChoice(){return timeFitChoice;}
    string* GetElementNames(){return elementNames;}
    //functions
    void CreateRequiredDataSets();
    void CreateShiftedBinEdges(Int_t numBins, Double_t binWidth, Double_t* binEdges);
    ~FitOption();
};

/// Creates any required data set such as bin number array, time fit end array, time fit start array, and bin width array.
/// These data sets are needed to be generated so that the things like the histograms can be generated.
void FitOption::CreateRequiredDataSets()
{
    binNumArr = new Int_t [numCycles];
    timeFitEndArr = new Double_t [numCycles];
    binWidthArr = new Double_t [numCycles];
    eventNumArr = new Double_t [numCycles];
    timeFitStartArr = new Double_t [numCycles];
    timeLengthArr = new Double_t [numCycles];
    fitStartBinArr = new Double_t [numCycles];
    fitEndBinArr = new Double_t [numCycles];
    binEdgesArr = new Double_t* [numCycles];

    //need to calculate bin width initally
    if(rebinChoice)
    {
        Double_t rebinSize = (Double_t) numBins;
        Double_t tempWidth;

        //bin data for Bateman histograms
        for(int i = 0; i < numCycles; i++)
        {
            timeFitEndArr[i] = timeRunEndSimulated;
            tempWidth = timeRunEndSimulated / rebinSize;
            binWidthArr[i] = tempWidth;
            binNumArr[i] = rebinSize;
            timeFitStartArr[i] = 0.0f;
            rebinSize = rebinSize + rebinBinInc;
            tempWidth = timeRunEndSimulated / rebinSize;
        }

        //bin data for integral histograms
        for(int i = 0; i < numCycles; i++)
        {
            tempWidth = binWidthArr[i];
            rebinSize = binNumArr[i];
            int intRebinSize = (int)rebinSize;
            binEdgesArr[i] = new Double_t[intRebinSize + 2];
            CreateShiftedBinEdges(rebinSize, tempWidth, binEdgesArr[i]);
        }
    }
    //need to calculate number of bins initally, change time fit based on number bins
    if(!rebinChoice && inputHistoExecutionType == 1)
    {
        Double_t fitEnd = timeRunEndSimulated;
        Double_t fitStart = 0.0f;
        Double_t initBinNum = fitEnd / binWidth;
        Double_t addedBins = 0.0f;
        Double_t totalBins = initBinNum;
        //data for Bateman histograms
        for(int i = 0; i < numCycles; i++)
        {
            timeFitEndArr[i] = fitEnd;
            binNumArr[i] = totalBins;
            binWidthArr[i] = binWidth;
            timeFitStartArr[i] = fitStart;

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

        //data for integral histograms
        for(int i = 0; i < numCycles; i++)
        {
            binEdgesArr[i] = new Double_t[binNumArr[i] + 2];
            CreateShiftedBinEdges(binNumArr[i], binWidthArr[i], binEdgesArr[i]);
        }
    }
    //for input histogram
    if(inputHistoExecutionType == 3 || inputHistoExecutionType == 2)
    {
        Double_t binWidth = inputHistoTimeEnd / inputHistoBinNum;
        Double_t fitStart = timeRunStartInput;
        Double_t fitEnd = timeRunEndInput;
        Double_t fitStartBinNum;
        Double_t fitEndBinNum;
        //data for input histogram
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
        //data for integral of input histogram
        for(int i = 0; i < numCycles; i++)
        {
            binEdgesArr[i] = new Double_t[binNumArr[i] + 2];
            CreateShiftedBinEdges(binNumArr[i], binWidth, binEdgesArr[i]);
        }
    }
    //Getting number of events generated for each different cycle
    if(eventNumChangeChoice)
    {
        cout << "Change factor " << eventChangeFactor << endl;
        eventNumArr[0] = events;
        for(int i = 1; i < numCycles; i++)
        {
            if(eventChangeType == 1)
            {
                eventNumArr[i] = TMath::Power(eventChangeFactor, i) * events;
            }else if(eventChangeType == 2)
            {
                eventNumArr[i] = (eventChangeFactor * i) + events;
            }
        }
    }
}

///Creates the bin edges arrays required for the creation of the integral histograms
void FitOption::CreateShiftedBinEdges(Int_t numBins, Double_t binWidth, Double_t* binEdges)
{
    Double_t tempBinEdge = 0.0;

    binEdges[0] = 0.0;
    binEdges[1] = (binWidth / 2.0);
    for(int i = 2; i < (numBins + 2); i++)
    {
        binEdges[i] = ((i-1) * binWidth) + (binWidth / 2.0);
    }
}

FitOption::~FitOption()
{
    delete [] binNumArr;
    delete [] binWidthArr;
    delete [] timeFitEndArr;
    delete [] timeFitStartArr;
    delete [] eventNumArr;
    delete [] timeLengthArr;
    delete [] fitStartBinArr;
    delete [] fitEndBinArr;
    for(int i = 0; i < numCycles; i++)
    {
        delete [] binEdgesArr[i];
    }
    delete [] binEdgesArr;
}

#endif