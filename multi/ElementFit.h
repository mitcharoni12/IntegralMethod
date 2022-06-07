#ifndef ELEMENTFIT_H
#define ELEMENTFIT_H

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "TF1.h"
#include "TH1.h"
#include "time.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"

#include "SingleElementFitValues.h"
#include "ChainFitValues.h"
#include "SingleCycleHistoHolder.h"
#include "CycleHistoHolder.h"
#include "SingleCycleCanvasHolder.h"
#include "CycleCanvasHolder.h"
#include "CycleIntegralHistoHolder.h"
#include "SingleCycleIntegralHistoHolder.h"
#include "ParameterValue.h"
#include "FitFunction.h"
#include "FitOption.h"
#include "Math/MinimizerOptions.h"

using namespace std;

typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);

/// This class can be considered the cornerstore for the entire program. This class is the class that generated the data then fits it.
/// After the generation and fitting, the fitted values are passed onto the Run class if multiple Runs/Cycles are chosen.
class ElementFit{
    public:
        ElementFit(Double_t (*batemanFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Double_t (**batemanFitFunctions)(Double_t*, Double_t*), Double_t (**integralFitFunctions)(Double_t*, Double_t*),
                   ParameterValue** paraVals, FitOption* fitOptions);
        ElementFit(Double_t (*batemanFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Double_t (**batemanFitFunctions)(Double_t*, Double_t*), Double_t (**integralFitFunctions)(Double_t*, Double_t*),
                   ParameterValue** paraVals, FitOption* fitOptions, TH1D* inputHistogram);
        ~ElementFit();
        //getter function
        Double_t GetElementParameters(int i){return paraVals[i]->GetDecayConst();}              ///< Returns decay constant values used to set starting point for fit functions
        FitOption* GetFitOptions(){return fitOptions;}                                          ///< Returns options the user chose to fit with
        ChainFitValues* GetBatemanFitValues(){return totalBatemanFitValues;}                ///< Returns fitted values for the fit of the total Bateman histogram
        ChainFitValues* GetIntegralFitValues(){return totalIntegralFitValues;}              ///< Returns fitted values for the fit of the total Integral histogram
        SingleElementFitValues* GetSingleBatemanFitValues(){return singleBatemanFitValues;}     ///< Returns fitted values for the fit of all the single Bateman histograms
        SingleElementFitValues* GetSingleIntegralFitValues(){return singleIntegralFitValues;}   ///< Returns fitted values for the fit of all the single Integral histograms
        //setter function
        void setNumRuns(Int_t numRuns){this->numRuns = numRuns;}
        void setNumCycles(Int_t numCycles){this->numCycles = numCycles;}

        void ChangeSeed();
        void CreateInputIntegralHisto();
        void CutInputHistos();
        void DisplayTotalFunctionParameters();
        void DisplayParameterLimits();
        void GenBatemanHistograms();
        //functions for single bateman data
        void CreateSingleBatemanFitFunctions(Int_t timeEnd);
        void CreateSingleBatemanHistoHolders();
        void DisplaySingleBatemanHistos(TCanvas** can);
        void DisplaySingleBatemanParameters();
        void DrawSingleBatemanIndividualHistos(SingleCycleCanvasHolder* singleBatemanCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex);
        void FitSingleBatemanHistos(Int_t cycleIndex, Int_t runIndex);
        void SetSingleBatemanFunctionParameters();
        void SetSingleBatemanParameterLimits();
        //functions for single integral data
        void CreateSingleIntegralFitFunctions(Int_t timeEnd);
        void GenSingleIntegralHisto();
        void CreateSingleIntegralHistoHolders();
        void DisplaySingleIntegralHisto(TCanvas** can);
        void DisplaySingleIntegralParameters();
        void DrawSingleIntegralIndividualHistos(SingleCycleCanvasHolder* singleIntegralCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex);
        void FitSingleIntegralHistos(Int_t cycleIndex, Int_t runIndex);
        void GenSingleIntegralHistos();
        void SetSingleIntegralFunctionParameters();
        void SetSingleIntegralParameterLimits();
        //functions for total bateman data
        void CreateTotalBatemanFitFunctions(Int_t timeEnd);
        void CreateTotalBatemanHistoHolders();
        void DisplayTotalBatemanHisto(TCanvas* can);
        void DisplayTotalBatemanParameters();
        void DrawTotalBatemanIndividualHistos(CycleCanvasHolder* batemanTotalCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex);
        void FitTotalBatemanHistos(Int_t cycleIndex, Int_t runIndex);
        void SetTotalBatemanFunctionParameters();
        void SetTotalBatemanParameterLimits();
        //functions for total integral data
        void CreateTotalIntegralFitFunctions(Int_t timeEnd);
        void GenTotalIntegralHisto();
        void CreateTotalIntegralHistoHolders();
        void DisplayTotalIntegralHisto(TCanvas* can);
        void DisplayTotalIntegralParameters();
        void DrawTotalIntegralIndividualHistos(CycleCanvasHolder* integralTotalCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex);
        void FitTotalIntegralHistos(Int_t cycleIndex, Int_t runIndex);
        void SetTotalIntegralFunctionParameters();
        void SetTotalIntegralParameterLimits();
    private:
        //private variables
        FitOption* fitOptions;                                                      ///< Contains fit options for program
        string* elementNames;                                                       ///< Contains element names for every element in the decay chain
        Int_t numEvents, numElements, numParameters, numRuns, numCycles, inputHistoExecutionType, singleElementDataChoice;
        bool multiSourceChoice;                                                     ///< True = generate multiple sets of histograms for the different cyles. False = generate a single set of histograms for the diiferent cyles.
        bool rebinChoice;                                                           ///< True = Have program execute with changing bin number between cylces and keeping time fit constant. False = Dont do rebin
        bool numEventChangeChoice;                                                  ///< True = Change number of events between cylces
        bool inputHistogramChoice;                                                  ///< True = dealing with input histogram.
        Int_t* binNumArr;                                                           ///< Contains an array of number of bins the histogram have between cycles.
        TF1* integralFunction, *batemanFunction;
        TF1** singleIntegralFitFunctions, **singleBatemanFitFunctions;
        CycleHistoHolder* batemanHisto;                                             ///< Contains every total Bateman histogram for the entire program.
        CycleIntegralHistoHolder* integralHisto;                                    ///< Contains every total integral histogram for the entire program.
        SingleCycleHistoHolder* singleBatemanHisto;                                 ///< Contains every single Bateman histogram for the entire program.
        SingleCycleIntegralHistoHolder* singleIntegralHisto;                        ///< Contains every single integral histogram for the entire program.
        TH1D* inputHistogram;                                                       ///< Contains the raw input histogram.
        //passed fit functions used to make the TF1 for fitting
        decayFunction* batemanFitFunctions, *integralFitFunctions;
        decayFunction passedBatemanFunction, passedIntegralFunction;
        TRandom3 rand;                                                              ///< Used for random number generation.
        Double_t timeRunStart;                                                      ///< Inital value for the start of fit time, can change between cycle.
        Double_t timeRunEnd;                                                        ///< Inital value for the end of fit time, can change between cycle.
        Double_t leaveOutStartBinsSim;                                              ///< How many bins to leave out of the start of the fit for the simulated histograms.
        Double_t leaveOutEndBinsSim;                                                ///< How many bins to leave out of the end of the fit for the simulated histograms.
        Double_t leaveOutStartBinsInput;                                            ///< How many bins to leave out of the start of the fit for the input histogram.
        Double_t leaveOutEndBinsInput;                                              ///< How many bins to leave out of the end of the fit for the input histogram.
        Double_t doubleNumEvents;                                                   ///< Number of events except its of data type Double_t.
        Double_t* randArr;                                                          ///< Array used for storing random values for event generation.
        Double_t* binWidth;                                                         ///< Holds bin widths between cycles.
        Double_t* timeEndArr;                                                       ///< Holds the time of fit end between cycles.
        Double_t* timeStartArr;                                                     ///< Holds the time of fit start between cycles.
        Double_t** binEdgesArr;                                                     ///< Holds the bin edges for the integral histograms.
        ChainFitValues* totalBatemanFitValues, *totalIntegralFitValues;
        SingleElementFitValues* singleBatemanFitValues, *singleIntegralFitValues;
        ParameterValue** paraVals;                                                  ///< Accepted values for all parameters in all the fit functions.
        Int_t globalSeedChanger = 0;                                                ///< Used for assigning in changing the seed.
        TCanvas* test = new TCanvas("test", "test", 500, 500);
        //TCanvas* test2 = new TCanvas("test2", "test2", 500, 500);
        clock_t Tclock;
};  

/// \brief The constructor for simualting data
///
/// Will take the fit fuctions for the program and the fit options for both the Run class and Cycle class. Dynamically allocates the fit parameter storages here.
/// Creates data required for creating the histograms between runs and cycles. Runs the call for creating all the histograms.
ElementFit::ElementFit(Double_t (*batemanFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Double_t (**batemanFitFunctions)(Double_t*, Double_t*), Double_t (**integralFitFunctions)(Double_t*, Double_t*),
                   ParameterValue** paraVals, FitOption* fitOptions)
{
    //setting variables
    this->fitOptions = fitOptions;
    this->numEvents = fitOptions->GetNumEvents();
    this->numRuns = fitOptions->GetNumRuns();
    this->numCycles = fitOptions->GetNumCycles();
    this->elementNames = fitOptions->GetElementNames();
    this->batemanFitFunctions = batemanFitFunctions;
    this->integralFitFunctions = integralFitFunctions;
    this->numElements = fitOptions->GetNumElements();
    this->numParameters = numElements*2;
    this->passedBatemanFunction = batemanFunc;
    this->passedIntegralFunction = integralFunc;
    this->timeRunEnd = fitOptions->GetTimeRunEndSimulated();
    this->paraVals = paraVals;
    this->multiSourceChoice = fitOptions->GetMultiSourceChoice();
    this->numEventChangeChoice = fitOptions->GetNumEventChangeChoice();
    this->rebinChoice = fitOptions->GetRebinChoice();
    this->leaveOutStartBinsSim = fitOptions->GetLeaveOutStartBinsSim();
    this->leaveOutEndBinsSim = fitOptions->GetLeaveOutEndBinsSim();
    this->leaveOutStartBinsInput = fitOptions->GetLeaveOutStartBinsInput();
    this->leaveOutEndBinsInput = fitOptions->GetLeaveOutEndBinsInput();
    this->inputHistoExecutionType = fitOptions->GetInputHistoExecutionType();
    this->singleElementDataChoice = fitOptions->GetSingleElementDataChoice();
    doubleNumEvents = (Double_t) numEvents;

    //creating fit parameter holders
    singleBatemanFitValues = new SingleElementFitValues(numElements);
    singleIntegralFitValues = new SingleElementFitValues(numElements);
    totalBatemanFitValues = new ChainFitValues(numElements);
    totalIntegralFitValues = new ChainFitValues(numElements);

    singleBatemanFitFunctions = new TF1* [numElements];
    singleIntegralFitFunctions = new TF1* [numElements];
    //setting values for randomization
    randArr = new Double_t [numElements];
    Tclock = clock();
    //creates required arrays for number of bins and time fit
    fitOptions->CreateRequiredDataSets();
    binWidth = fitOptions->GetBinWidthArr();
    timeEndArr = fitOptions->GetTimeFitEndArr();
    timeStartArr = fitOptions->GetTimeFitStartArr();
    binNumArr = fitOptions->GetBinNumArr();
    binEdgesArr = fitOptions->GetBinEdges();
    //generates all the histo objects and fills each one with data
    CreateTotalBatemanHistoHolders();
    CreateTotalIntegralHistoHolders();
    if(singleElementDataChoice == 2)
    {
        CreateSingleBatemanHistoHolders();
        CreateSingleIntegralHistoHolders();
    }
    GenBatemanHistograms();
    GenTotalIntegralHisto();
    if(singleElementDataChoice == 2)
    {
        GenSingleIntegralHistos();
    }
    //removes restrictions on fitting function calls and itterations
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000000);
}


/// \brief Constructor for input histogram
ElementFit::ElementFit(Double_t (*batemanFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Double_t (**batemanFitFunctions)(Double_t*, Double_t*), Double_t (**integralFitFunctions)(Double_t*, Double_t*),
                   ParameterValue** paraVals, FitOption* fitOptions, TH1D* inputHistogram)
{
    //setting variables
    this->fitOptions = fitOptions;
    this->numRuns = fitOptions->GetNumRuns();
    this->numCycles = fitOptions->GetNumCycles();
    this->elementNames = fitOptions->GetElementNames();
    this->batemanFitFunctions = batemanFitFunctions;
    this->integralFitFunctions = integralFitFunctions;
    this->numElements = fitOptions->GetNumElements();
    this->numParameters = numElements*2;
    this->passedBatemanFunction = batemanFunc;
    this->passedIntegralFunction = integralFunc;
    this->timeRunEnd = fitOptions->GetTimeRunEndInput();
    this->paraVals = paraVals;
    this->multiSourceChoice = fitOptions->GetMultiSourceChoice();
    this->rebinChoice = fitOptions->GetRebinChoice();
    this->leaveOutStartBinsSim = fitOptions->GetLeaveOutStartBinsSim();
    this->leaveOutEndBinsSim = fitOptions->GetLeaveOutEndBinsSim();
    this->leaveOutStartBinsInput = fitOptions->GetLeaveOutStartBinsInput();
    this->leaveOutEndBinsInput = fitOptions->GetLeaveOutEndBinsInput();
    this->inputHistogram = inputHistogram;
    this->inputHistoExecutionType = fitOptions->GetInputHistoExecutionType();
    this->singleElementDataChoice = fitOptions->GetSingleElementDataChoice();
    doubleNumEvents = (Double_t) numEvents;

    //creating fit parameter holders
    singleBatemanFitValues = new SingleElementFitValues(numElements);
    singleIntegralFitValues = new SingleElementFitValues(numElements);
    totalBatemanFitValues = new ChainFitValues(numElements);
    totalIntegralFitValues = new ChainFitValues(numElements);

    singleBatemanFitFunctions = new TF1* [numElements];
    singleIntegralFitFunctions = new TF1* [numElements];
    //setting values for randomization
    randArr = new Double_t [numElements];
    Tclock = clock();
    //creates required arrays for number of bins and time fit
    fitOptions->CreateRequiredDataSets();
    binWidth = fitOptions->GetBinWidthArr();
    timeEndArr = fitOptions->GetTimeFitEndArr();
    timeStartArr = fitOptions->GetTimeFitStartArr();
    binNumArr = fitOptions->GetBinNumArr();
    binEdgesArr = fitOptions->GetBinEdges();
    CreateTotalBatemanHistoHolders();
    CutInputHistos();
    CreateInputIntegralHisto();
    //removes restrictions on fitting function calls and itterations
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000000);
}

ElementFit::~ElementFit()
{
    if(inputHistoExecutionType == 1)
    {    
        delete totalBatemanFitValues;
        delete totalIntegralFitValues;
        delete [] randArr;
        delete batemanHisto;
        delete integralHisto;
        if(singleElementDataChoice == 2)
        {
            delete singleBatemanFitValues;
            delete singleIntegralFitValues;
            delete [] singleBatemanFitFunctions;
            delete [] singleIntegralFitFunctions;
            delete singleBatemanHisto;
            delete singleIntegralHisto;
        }
    }else if(inputHistoExecutionType == 2 || inputHistoExecutionType == 3)
    {
        delete totalBatemanFitValues;
        delete totalIntegralFitValues;
        delete [] randArr;
        if(singleElementDataChoice == 2)
        {
            delete singleBatemanFitValues;
            delete singleIntegralFitValues;
            delete [] singleBatemanFitFunctions;
            delete [] singleIntegralFitFunctions;
        }
    }
}

/// \brief changes seed for random number generator
void ElementFit::ChangeSeed()
{
    Double_t seederSeed;
    Double_t randNumber;

    randNumber = rand.Uniform();
    seederSeed = Tclock + globalSeedChanger + randNumber;

    rand.SetSeed(Tclock + globalSeedChanger + seederSeed);
    globalSeedChanger++;
}

/// \brief Creates the integral histos for the input histogram.
void ElementFit::CreateInputIntegralHisto()
{
    string histoName = "Total Input Integral Histo";
    TH1D* tempBatemanHisto;
    TH1D* tempIntegralHisto;

    //creates the integral histogram objects
    integralHisto = new CycleIntegralHistoHolder(numCycles, numRuns, histoName, binNumArr, binEdgesArr);

    //fills integral histo
    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        tempBatemanHisto = batemanHisto->GetAHisto(cycleIndex, 0);
        tempIntegralHisto = integralHisto->GetAHisto(cycleIndex, 0);

        tempIntegralHisto->SetBinContent(1, 0.0);
        tempIntegralHisto->SetBinContent(2, tempBatemanHisto->GetBinContent(2));
        for(int binIndex = 3; binIndex < (binNumArr[cycleIndex] + 1); binIndex++)
        {
            tempIntegralHisto->SetBinContent(binIndex, tempIntegralHisto->GetBinContent(binIndex-1) + tempBatemanHisto->GetBinContent(binIndex-1));
        }
    }
}

/// \brief Cuts the original input histogram bins and puts them in a new histogram.
void ElementFit::CutInputHistos()
{
    TH1D* tempHisto;
    Double_t* tempBinStartArr = fitOptions->GetFitStartBinArr();
    Double_t* tempBinEndArr = fitOptions->GetFitEndBinArr();
    Double_t tempStartBin, tempBinEnd, tempBinValue;
    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        tempHisto = batemanHisto->GetAHisto(cycleIndex, 0);
        tempStartBin = tempBinStartArr[cycleIndex];
        for(int binIndex = 1; binIndex <= binNumArr[cycleIndex]; binIndex++)
        {
            tempBinValue = inputHistogram->GetBinContent(tempStartBin + binIndex - 2);
            tempHisto->SetBinContent(binIndex, tempBinValue);
        }
    }
}


/// \brief Creates all single bateman histograms for the entre program.
void ElementFit::CreateSingleBatemanHistoHolders()
{
    string histoName;
    histoName = "Single Bateman Histo";
    singleBatemanHisto = new SingleCycleHistoHolder(numCycles, numElements, numRuns, histoName, binNumArr, timeEndArr, elementNames);
}

/// \brief Creates all single integral histograms for the entre program.
void ElementFit::CreateSingleIntegralHistoHolders()
{
    string histoName;
    histoName = "Single Integral Histo";
    singleIntegralHisto = new SingleCycleIntegralHistoHolder(numCycles, numElements, numRuns, histoName, binNumArr, binEdgesArr, elementNames);
}

/// \brief Creates all total bateman histograms for the entre program.
void ElementFit::CreateTotalBatemanHistoHolders()
{
    string histoName;
    if(inputHistoExecutionType == 1)
    {
        histoName = "Total Bateman Histo";
        batemanHisto = new CycleHistoHolder(numCycles, numRuns, histoName, binNumArr, timeEndArr);
    }else if(inputHistoExecutionType == 2 || inputHistoExecutionType == 3)
    {
        Double_t* timeLengthArr = fitOptions->GetTimeLengthArr();
        histoName = "Total Input Bateman Histo";
        batemanHisto = new CycleHistoHolder(numCycles, numRuns, histoName, binNumArr, timeLengthArr);
    }
}

/// \brief Creates all total integral histograms for the entre program.
void ElementFit::CreateTotalIntegralHistoHolders()
{
    string histoName;
    histoName = "Total Integral Histo";
    integralHisto = new CycleIntegralHistoHolder(numCycles, numRuns, histoName, binNumArr, binEdgesArr);
}

/// \brief creates the TF1 objects for the single Bateman fit functions 
void ElementFit::CreateSingleBatemanFitFunctions(Int_t timeEnd)
{
    for(int i = 0; i < numElements; i++)
    {
        singleBatemanFitFunctions[i] = new TF1((elementNames[i] + "RegSingFunc").c_str(), batemanFitFunctions[i], 0., timeEnd, (i+1)*2);
    }
}

/// \brief creates the TF1 objects for the single integral fit functions 
void ElementFit::CreateSingleIntegralFitFunctions(Int_t timeEnd)
{
    for(int i = 0; i < numElements; i++)
    {
        singleIntegralFitFunctions[i] = new TF1((elementNames[i] + "InteSingFunc").c_str(), integralFitFunctions[i], 0., timeEnd, (i+1)*2);
    }
}

/// \brief creates the TF1 objects for the total Bateman fit function
void ElementFit::CreateTotalBatemanFitFunctions(Int_t timeEnd)
{
    batemanFunction = new TF1("TotalBatemanFunction", passedBatemanFunction, 0., timeEnd, numParameters);
}

/// \brief creates the TF1 objects for the total integral fit function
void ElementFit::CreateTotalIntegralFitFunctions(Int_t timeEnd)
{
    integralFunction = new TF1("TotalIntegralFunction", passedIntegralFunction, 0., timeEnd, numParameters);
}

/// \brief Displays fitted parameters of the single bateman for the single run option for the program
void ElementFit::DisplaySingleBatemanParameters()
{
    cout << "BATEMAN SINGLE FIT PARAMETERS/ERRORS" << endl << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << " Single Fit Values" << endl;
        for(int k = 0; k < i+1; k++)
        {
            cout << elementNames[k] << ": \tHalf Life: " << singleBatemanFitValues->GetAnHalfLife(i, k) << "s" << endl << 
            "\tHalf Life Error: " << singleBatemanFitValues->GetAnHalfLifeError(i,k) << "s" << endl;
            cout << "\tN0: " << singleBatemanFitValues->GetAnN0(i, k) << "" << endl << 
            "\tN0 Error: " << singleBatemanFitValues->GetAnN0Error(i,k) << "" << endl;
            cout << endl;
        }
    }
}

/// \brief Displays fitted parameters of the single integral for the single run option for the program
void ElementFit::DisplaySingleIntegralParameters()
{
    cout << endl << "INTEGRAL SINGLE FIT PARAMETERS/ERRORS" << endl << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << " Single Fit Values" << endl;
        for(int k = 0; k < i+1; k++)
        {
            cout << elementNames[k] << ": \tHalf Life: " << singleIntegralFitValues->GetAnHalfLife(i, k) << "s" << endl << 
            "\tHalf Life Error: " << singleIntegralFitValues->GetAnHalfLifeError(i,k) << "s" << endl;
            cout << "\tN0: " << singleIntegralFitValues->GetAnN0(i, k) << "" << endl << 
            "\tN0 Error: " << singleIntegralFitValues->GetAnN0Error(i,k) << "" << endl;
            cout << endl;
        }
    }
}

/// \brief Displays fitted parameters of the total bateman for the single run option for the program
void ElementFit::DisplayTotalBatemanParameters()
{  
    cout << endl << "BATEMAN FIT PARAMETERS/ERRORS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ": " << "\tHalfLife: " << totalBatemanFitValues->GetAnHalfLife(i) << "s" << endl << 
        "\tError: " << totalBatemanFitValues->GetAnHalfLifeError(i) << "s" << endl;
        cout << "\tN0: " << totalBatemanFitValues->GetAnN0(i) << "" << endl << 
        "\tError: " << totalBatemanFitValues->GetAnN0Error(i) << "" << endl << endl;
    }
}

/// \brief Displays fitted parameters of the total integral for the single run option for the program
void ElementFit::DisplayTotalIntegralParameters()
{
    cout << "INTEGRAL FIT PARAMETERS/ERRORS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ": " << "\tHalfLife: " << totalIntegralFitValues->GetAnHalfLife(i) << "s" << endl << 
        "\tError: " << totalIntegralFitValues->GetAnHalfLifeError(i) << "s" << endl;
        cout << "\tN0: " << totalIntegralFitValues->GetAnN0(i) << "" << endl << 
        "\tError: " << totalIntegralFitValues->GetAnN0Error(i) << "" << endl << endl;
    }
}

/// \brief Displays single Bateman histograms for the single run option for the program.
void ElementFit::DisplaySingleBatemanHistos(TCanvas** canvas)
{
    for(int i = 0; i < numElements; i++)
    {
        canvas[i]->cd();
        singleBatemanHisto->GetAHisto(0, 0, i)->Draw();
    }
}

/// \brief Displays single integral histograms for the single run option for the program.
void ElementFit::DisplaySingleIntegralHisto(TCanvas** canvas)
{
    for(int i = 0; i < numElements; i++)
    {
        canvas[i]->cd();
        singleIntegralHisto->GetAHisto(0, 0, i)->Draw();
    }
}

/// \brief Displays total Bateman histogram for the single run option for the program
void ElementFit::DisplayTotalBatemanHisto(TCanvas* canvas)
{
    canvas->cd();
    batemanHisto->GetAHisto(0, 0)->Draw();
}

/// \brief Displays total Bateman histogram for the single run option for the program
void ElementFit::DisplayTotalIntegralHisto(TCanvas* canvas)
{
    canvas->cd();
    integralHisto->GetAHisto(0, 0)->Draw();
}

/// \brief Displays the single Bateman individual histogram for any execution type of the program.
void ElementFit::DrawSingleBatemanIndividualHistos(SingleCycleCanvasHolder* singleBatemanCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex)
{
    upperRunIndex = upperRunIndex + 1;
    upperCycleIndex = upperCycleIndex + 1;

    for(int cycleIndex = lowerCycleIndex; cycleIndex < upperCycleIndex; cycleIndex++)
    {
        for(int runIndex = lowerRunIndex; runIndex < upperRunIndex; runIndex++)
        {
            for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
            {
                singleBatemanCanvases->GetACanvas(cycleIndex, runIndex, elementIndex)->cd();
                singleBatemanHisto->GetAHisto(cycleIndex, runIndex, elementIndex)->Draw();
            }
        }
    }
}

/// \brief Displays the single Bateman individual histogram for any execution type of the program.
void ElementFit::DrawSingleIntegralIndividualHistos(SingleCycleCanvasHolder* singleIntegralCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex)
{
    upperRunIndex = upperRunIndex + 1;
    upperCycleIndex = upperCycleIndex + 1;

    for(int cycleIndex = lowerCycleIndex; cycleIndex < upperCycleIndex; cycleIndex++)
    {
        for(int runIndex = lowerRunIndex; runIndex < upperRunIndex; runIndex++)
        {
            for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
            {
                singleIntegralCanvases->GetACanvas(cycleIndex, runIndex, elementIndex)->cd();
                singleIntegralHisto->GetAHisto(cycleIndex, runIndex, elementIndex)->Draw();
            }
        }
    }
}

/// \brief Displays the total Bateman individual histogram for any execution type of the program
void ElementFit::DrawTotalBatemanIndividualHistos(CycleCanvasHolder* batemanTotalCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex)
{
    upperRunIndex = upperRunIndex + 1;
    upperCycleIndex = upperCycleIndex + 1;

    for(int cycleIndex = lowerCycleIndex; cycleIndex < upperCycleIndex; cycleIndex++)
    {
        for(int runIndex = lowerRunIndex; runIndex < upperRunIndex; runIndex++)
        {
            cout << "Getting cycle: " << cycleIndex << " Run: " << runIndex << endl;
            batemanTotalCanvases->GetACanvas(cycleIndex, runIndex)->cd();
            batemanHisto->GetAHisto(cycleIndex, runIndex)->Draw();
        }
    }
}

/// \brief Displays the total Bateman individual histogram for any execution type of the program
void ElementFit::DrawTotalIntegralIndividualHistos(CycleCanvasHolder* integralTotalCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex)
{
    upperRunIndex = upperRunIndex + 1;
    upperCycleIndex = upperCycleIndex + 1;

    for(int cycleIndex = lowerCycleIndex; cycleIndex < upperCycleIndex; cycleIndex++)
    {
        for(int runIndex = lowerRunIndex; runIndex < upperRunIndex; runIndex++)
        {
            integralTotalCanvases->GetACanvas(cycleIndex, runIndex)->cd();
            integralHisto->GetAHisto(cycleIndex, runIndex)->Draw();
        }
    }
}

/// \brief Fits single Bateman histograms at specified cycle and run index
void ElementFit::FitSingleBatemanHistos(Int_t cycleIndex, Int_t runIndex)
{
    Double_t valueN0, errorN0, valueDecayConst, errorDecayConst, valueHalfLife, errorHalfLife, startFitOffset, startFit, endFitOffset, endFit;
    TFitResultPtr fitResult = nullptr;
    TH1D* tempSingleBatemanHisto;

    timeRunStart = timeStartArr[cycleIndex];
    timeRunEnd = timeEndArr[cycleIndex];
    CreateSingleBatemanFitFunctions(timeEndArr[cycleIndex]);
    SetSingleBatemanFunctionParameters();
    SetSingleBatemanParameterLimits();

    startFitOffset = leaveOutStartBinsSim * binWidth[cycleIndex];
    startFit = startFitOffset + timeRunStart;
    endFitOffset = leaveOutEndBinsSim * binWidth[cycleIndex];
    endFit = timeRunEnd - endFitOffset;

    //loop for fitting and storing value for elements
    for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
    {
        tempSingleBatemanHisto = singleBatemanHisto->GetAHisto(cycleIndex, runIndex, elementIndex);

        cout << "FITTING SINGLE BATEMAN " << elementNames[elementIndex] << " CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;
        fitResult = tempSingleBatemanHisto->Fit(singleBatemanFitFunctions[elementIndex], "QLSMULTITHREAD", "", startFit, endFit);

        for(int subElementIndex = 0; subElementIndex < (elementIndex+1); subElementIndex++)
        {
            //storing values for bateman part of function
            valueN0 = singleBatemanFitFunctions[elementIndex]->GetParameter((subElementIndex*2));
            errorN0 = singleBatemanFitFunctions[elementIndex]->GetParError((subElementIndex*2));
            valueDecayConst = singleBatemanFitFunctions[elementIndex]->GetParameter((subElementIndex*2)+1);
            errorDecayConst = singleBatemanFitFunctions[elementIndex]->GetParError((subElementIndex*2)+1);
            valueHalfLife = log(2)/(valueDecayConst);
            errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));
            //if fait fails
            if(fitResult->Status() != 0)
            {
                valueHalfLife = 0.000000000000;
            }

            singleBatemanFitValues->SetAnN0(elementIndex, subElementIndex, valueN0);
            singleBatemanFitValues->SetAnN0Error(elementIndex, subElementIndex, errorN0);
            singleBatemanFitValues->SetAnHalfLife(elementIndex, subElementIndex, valueHalfLife);
            singleBatemanFitValues->SetAnHalfLifeError(elementIndex, subElementIndex, errorHalfLife);
        }
    }
}

/// \brief Fits single integral histos at specified cycle and run index
void ElementFit::FitSingleIntegralHistos(Int_t cycleIndex, Int_t runIndex)
{
    Double_t valueN0, errorN0, valueDecayConst, errorDecayConst, valueHalfLife, errorHalfLife, startFitOffset, startFit, endFitOffset, endFit;
    TH1D* tempSingleBatemanHisto, *tempSingleHisto;
    TFitResultPtr fitResult = nullptr;

    timeRunStart = timeStartArr[cycleIndex];
    timeRunEnd = timeEndArr[cycleIndex];
    CreateSingleIntegralFitFunctions(timeEndArr[cycleIndex]);
    SetSingleIntegralFunctionParameters();
    SetSingleIntegralParameterLimits();

    startFitOffset = leaveOutStartBinsSim * binWidth[cycleIndex];
    startFit = startFitOffset + timeRunStart;
    endFitOffset = leaveOutEndBinsSim * binWidth[cycleIndex];
    endFit = timeRunEnd - endFitOffset;
    startFit = startFit + (binWidth[cycleIndex]/2.0);
    endFit = endFit + (binWidth[cycleIndex]/2.0);

    //loop for fitting and storing value for elements
    for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
    {
        tempSingleHisto = singleIntegralHisto->GetAHisto(cycleIndex, runIndex, elementIndex);

        cout << "FITTING SINGLE INTEGRAL " << elementNames[elementIndex] << " CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;
        fitResult = tempSingleHisto->Fit(singleIntegralFitFunctions[elementIndex], "QLSMULTITHREAD", "", startFit, endFit);

        for(int elementSubIndex = 0; elementSubIndex < (elementIndex+1); elementSubIndex++)
        {
            //storing values for integral part of function
            valueN0 = singleIntegralFitFunctions[elementIndex]->GetParameter((elementSubIndex*2));
            errorN0 = singleIntegralFitFunctions[elementIndex]->GetParError((elementSubIndex*2));
            valueDecayConst = singleIntegralFitFunctions[elementIndex]->GetParameter((elementSubIndex*2)+1);
            errorDecayConst = singleIntegralFitFunctions[elementIndex]->GetParError((elementSubIndex*2)+1);
            valueHalfLife = log(2)/(valueDecayConst);
            errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));
            //if fit fails
            if(fitResult->Status() != 0)
            {
                valueHalfLife = 0.000000000000;
            }

            singleIntegralFitValues->SetAnN0(elementIndex, elementSubIndex, valueN0);
            singleIntegralFitValues->SetAnN0Error(elementIndex, elementSubIndex, errorN0);
            singleIntegralFitValues->SetAnHalfLife(elementIndex, elementSubIndex, valueHalfLife);
            singleIntegralFitValues->SetAnHalfLifeError(elementIndex, elementSubIndex, errorHalfLife);
        }
    }
}

/// \brief Fits total Bateman histogram at specified cycle and run index
void ElementFit::FitTotalBatemanHistos(Int_t cycleIndex, Int_t runIndex)
{
    //dynamic array
    Double_t valueN0, errorN0, valueDecayConst, errorDecayConst, valueHalfLife, errorHalfLife, startFitOffset, startFit, endFitOffset, endFit;
    TFitResultPtr fitResult = nullptr;
    TH1D* tempHisto;

    timeRunStart = timeStartArr[cycleIndex];
    if(inputHistoExecutionType == 1)
    {
        timeRunEnd = timeEndArr[cycleIndex];
    }else if(inputHistoExecutionType == 2 || inputHistoExecutionType == 3)
    {
        timeRunEnd = fitOptions->GetTimeLengthArr()[cycleIndex];
    }
    CreateTotalBatemanFitFunctions(timeEndArr[cycleIndex]);
    SetTotalBatemanFunctionParameters();
    SetTotalBatemanParameterLimits();

    startFitOffset = leaveOutStartBinsSim * binWidth[cycleIndex];
    startFit = startFitOffset + timeRunStart;
    endFitOffset = leaveOutEndBinsSim * binWidth[cycleIndex];
    endFit = timeRunEnd - endFitOffset;

    tempHisto = batemanHisto->GetAHisto(cycleIndex, runIndex);

    cout << "FITTING TOTAL BATEMAN CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;
    //retry fit if it fails
    for(int i = 0; i < 1; i++)
    {
        fitResult = tempHisto->Fit(batemanFunction, "LSMULTITHREAD", "", startFit, endFit);
    }
    
    for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
    {   
        valueN0 = batemanFunction->GetParameter((elementIndex*2));
        errorN0 = batemanFunction->GetParError((elementIndex*2));
        valueDecayConst = batemanFunction->GetParameter((elementIndex*2)+1);
        errorDecayConst = batemanFunction->GetParError((elementIndex*2)+1);
        valueHalfLife = log(2)/(valueDecayConst);
        errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));
        if((inputHistoExecutionType == 2|| inputHistoExecutionType == 3) && (paraVals[elementIndex]->GetFixDecayConst() == false))
        {
            valueDecayConst = valueDecayConst / 1e-8;
            errorDecayConst = errorDecayConst / 1e-8;
            valueHalfLife = log(2)/(valueDecayConst);
            errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));
        }
        //if fit fails
        if(fitResult->Status() != 0)
        {
            valueHalfLife = 0.000000000000;
        }

        //used to get either the N0 value or half life value
        totalBatemanFitValues->SetAnN0(elementIndex, valueN0);
        totalBatemanFitValues->SetAnN0Error(elementIndex, errorN0);
        totalBatemanFitValues->SetAnHalfLife(elementIndex, valueHalfLife);
        totalBatemanFitValues->SetAnHalfLifeError(elementIndex, errorHalfLife);
    }
}

/// \brief Fits total integral histogram at specified cycle and run index
void ElementFit::FitTotalIntegralHistos(Int_t cycleIndex, Int_t runIndex)
{
    Double_t valueN0, errorN0, valueDecayConst, errorDecayConst, valueHalfLife, errorHalfLife, startFitOffset, startFit, endFitOffset, endFit;
    TFitResultPtr fitResult = nullptr;
    TH1D* tempHisto;

    timeRunStart = timeStartArr[cycleIndex];
    if(inputHistoExecutionType == 1)
    {
        timeRunEnd = timeEndArr[cycleIndex];
    }else if(inputHistoExecutionType == 2 || inputHistoExecutionType == 3)
    {
        timeRunEnd = fitOptions->GetTimeLengthArr()[cycleIndex];
    }
    CreateTotalIntegralFitFunctions(timeRunEnd);
    SetTotalIntegralFunctionParameters();
    SetTotalIntegralParameterLimits();

    startFitOffset = leaveOutStartBinsSim * binWidth[cycleIndex];
    startFit = startFitOffset + timeRunStart;
    endFitOffset = leaveOutEndBinsSim * binWidth[cycleIndex];
    endFit = timeRunEnd - endFitOffset;
    startFit = startFit + (binWidth[cycleIndex]/2.0);
    endFit = endFit + (binWidth[cycleIndex]/2.0);

    tempHisto = integralHisto->GetAHisto(cycleIndex, runIndex);

    cout << "FITTING TOTAL INTEGRAL CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;
    fitResult = tempHisto->Fit(integralFunction, "LSMULTITHREAD", "", startFit, endFit);

    for(int i = 0; i < numElements; i++)
    {   
        valueN0 = integralFunction->GetParameter((i*2));
        errorN0 = integralFunction->GetParError((i*2));
        valueDecayConst = integralFunction->GetParameter((i*2)+1);
        errorDecayConst = integralFunction->GetParError((i*2)+1);
        valueHalfLife = log(2)/(valueDecayConst);
        errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));
        if((inputHistoExecutionType == 2 || inputHistoExecutionType == 3) && (paraVals[i]->GetFixDecayConst() == false))
        {
            valueDecayConst = valueDecayConst / 1e-8;
            errorDecayConst = errorDecayConst / 1e-8;
            valueHalfLife = log(2)/(valueDecayConst);
            errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));
        }
        //if fit fails
        if(fitResult->Status() != 0)
        {
            valueHalfLife = 0.000000000000;
        }

        //used to get either the N0 value or half life value
        totalIntegralFitValues->SetAnN0(i, valueN0);
        totalIntegralFitValues->SetAnN0Error(i, errorN0);
        totalIntegralFitValues->SetAnHalfLife(i, valueHalfLife);
        totalIntegralFitValues->SetAnHalfLifeError(i, errorHalfLife);
    }
}

/// \brief Generates all the single integral histograms from the single Bateman histograms.
///
/// Generates all the single integral histograms from the single Bateman histograms.
/// Each individual bin(i) in the integral histogram is the sum of the integral bin i-1 plus the bin i from the regular histogram.
void ElementFit::GenSingleIntegralHistos()
{
    TH1D* tempIntegralHisto, *tempBatemanHisto;
    Int_t tempNumBins;

    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        tempNumBins = binNumArr[cycleIndex];
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
            {   
                tempIntegralHisto = singleIntegralHisto->GetAHisto(cycleIndex, runIndex, elementIndex);
                tempBatemanHisto = singleBatemanHisto->GetAHisto(cycleIndex, runIndex, elementIndex);

                tempIntegralHisto->SetBinContent(1, 0.0);
                tempIntegralHisto->SetBinContent(2, tempBatemanHisto->GetBinContent(1));
                for(int binIndex = 3; binIndex <= (tempNumBins+1); binIndex++)
                {
                    tempIntegralHisto->SetBinContent(binIndex, tempIntegralHisto->GetBinContent(binIndex-1) + tempBatemanHisto->GetBinContent(binIndex-1));
                }
            }
        }
    }
}

/// \brief Generates the total and single Bateman histograms
///
/// Generates the total and single Bateman histograms. Has different methods of generation for different program execution types.
/// Still not sure why generation works. This is the single most important function in the entire program.
void ElementFit::GenBatemanHistograms()
{
    TH1D* tempHisto;
    TH1D** singleTempHisto;
    singleTempHisto = new TH1D* [numElements];
    Double_t hold = 0.0f;
    Double_t stack = 0.0f;
    //case for generating the single histogram for the single source histogram choice
    if(!rebinChoice && !multiSourceChoice)
    {
        TH1D* tempHisto;
        TH1D** tempSingleHisto;
        tempSingleHisto = new TH1D* [numElements];

        ChangeSeed();

        //generates events for the single and total histograms symotaniously
        for(int i = 0; i < numEvents; i++)
        {
            //generating the times the events occured
            for(int j = 0; j < numElements; j++)
            {
                hold = rand.Uniform();
                randArr[j] = (-TMath::Log(hold)) / (paraVals[j]->GetDecayConst());
            }
            
            //putting events in repsective histograms
            for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
            {
                stack += randArr[elementIndex];
                for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
                {
                    tempHisto = batemanHisto->GetAHisto(cycleIndex, 0);
                    if(singleElementDataChoice == 2)
                    {
                        singleTempHisto[elementIndex] = singleBatemanHisto->GetAHisto(cycleIndex, 0, elementIndex);
                        singleTempHisto[elementIndex]->Fill(stack);
                    }
                    tempHisto->Fill(stack);
                }
            }
            stack = 0.0f; 
        }

        delete [] tempSingleHisto;
    //case for multiple histogram for the multiple source histogram choice.
    }else if(multiSourceChoice)
    {
        //data generated in run index 0 of cycle 0 must be identical to run index 0 of cycle 1 to see effects of rebinning and time fit change.
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            //change seed for generation between each run
            ChangeSeed();
            //generates events for the single and total histograms symotaniously
            for(int i = 0; i < numEvents; i++)
            {
                //generating the times the events occured
                for(int j = 0; j < numElements; j++)
                {
                    hold = rand.Uniform();
                    randArr[j] = (-TMath::Log(hold)) / (paraVals[j]->GetDecayConst());
                }
                //putting events in repsective histograms
                for(int k = 0; k < numElements; k++)
                {
                    stack += randArr[k];
                    //puts same events in all the cycles
                    //Ex: if we have 10 cycles, all the run 0's of all the cycles will have the same exact events.
                    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
                    {
                        tempHisto = batemanHisto->GetAHisto(cycleIndex, runIndex);
                        if(singleElementDataChoice == 2)
                        {
                            singleTempHisto[k] = singleBatemanHisto->GetAHisto(cycleIndex, runIndex, k);
                            singleTempHisto[k]->Fill(stack);
                        }
                        tempHisto->Fill(stack);
                    }
                }
                stack = 0.0f;
            }
        }
    //case for changing event number between cycles
    }else if(numEventChangeChoice)
    {
        Double_t* eventNums = fitOptions->GetEventNumArr();
    }

    delete [] singleTempHisto;
}

/// \brief Generates all total integral histograms for total Bateman histograms
///
/// Generates all total integral histograms for total Bateman histograms.
/// Each individual bin(i) in the integral histogram is the sum of the integral bin i-1 plus the bin i from the regular histogram.
void ElementFit::GenTotalIntegralHisto()
{
    //resets histogram so it can be used between runs
    TH1D* tempIntegralHisto, *tempBatemanHisto;
    Int_t tempNumBins;

    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        tempNumBins = binNumArr[cycleIndex];
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            tempIntegralHisto = integralHisto->GetAHisto(cycleIndex, runIndex);
            tempBatemanHisto = batemanHisto->GetAHisto(cycleIndex, runIndex);

            tempIntegralHisto->SetBinContent(1, 0.0);
            tempIntegralHisto->SetBinContent(2, tempBatemanHisto->GetBinContent(1));
            for(int binIndex = 3; binIndex <= (tempNumBins+1); binIndex++)
            {
                tempIntegralHisto->SetBinContent(binIndex, tempBatemanHisto->GetBinContent(binIndex-1) + tempIntegralHisto->GetBinContent(binIndex-1));
            }
        }
    }
}

/// \brief Sets the parameters for the single Bateman fit functions
void ElementFit::SetSingleBatemanFunctionParameters()
{
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k <= i; k++)
        {
            singleBatemanFitFunctions[i]->SetParameter((k*2), paraVals[k]->GetBatemanN0());
            singleBatemanFitFunctions[i]->SetParameter(((k*2)+1), paraVals[k]->GetDecayConst());
        }
    }
}

/// \brief Sets the parameters for the single Bateman fit functions
void ElementFit::SetSingleIntegralFunctionParameters()
{
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k <= i; k++)
        {
            singleIntegralFitFunctions[i]->SetParameter((k*2), paraVals[k]->GetIntegralN0());
            singleIntegralFitFunctions[i]->SetParameter(((k*2)+1), paraVals[k]->GetDecayConst());
        }
    }
}

/// \brief sets function parameters for the total bateman function
void ElementFit::SetTotalBatemanFunctionParameters()
{
    for(int i = 0; i < numElements; i++)
    {
        //have to have the if statments to account for the fact that the value could be a range
        if(inputHistoExecutionType == 1)
        {
            batemanFunction->SetParameter((i*2), paraVals[i]->GetBatemanN0());
            batemanFunction->SetParameter((i*2)+1, paraVals[i]->GetDecayConst());
        }else if(inputHistoExecutionType == 2 || inputHistoExecutionType == 3)
        {
            batemanFunction->SetParameter((i*2), paraVals[i]->GetBatemanN0());
            batemanFunction->SetParameter((i*2)+1, paraVals[i]->GetDecayConst10Ns());
        }
    }
}

/// \brief sets function parameters for total integral function
void ElementFit::SetTotalIntegralFunctionParameters()
{
    for(int i = 0; i < numElements; i++)
    {
        //have to have the if statments to account for the fact that the value could be a range
        if(inputHistoExecutionType == 1)
        {
            integralFunction->SetParameter((i*2), paraVals[i]->GetIntegralN0());
            integralFunction->SetParameter((i*2)+1, paraVals[i]->GetDecayConst());
        }else if(inputHistoExecutionType == 2 || inputHistoExecutionType == 3)
        {
            integralFunction->SetParameter((i*2), paraVals[i]->GetIntegralN0());
            integralFunction->SetParameter((i*2)+1, paraVals[i]->GetDecayConst10Ns());
        }
    }
}

/// \brief Displays the parameters used to set the fit functions
void ElementFit::DisplayTotalFunctionParameters()
{
    cout << "INITIAL FIT FUNCTION PARAMETERS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ":\tBateman init value: " << paraVals[i]->GetBatemanN0();
        cout << elementNames[i] << "\tIntegral init value: " << paraVals[i]->GetIntegralN0();
        cout << "\tHalf Life: " << paraVals[i]->GetHalfLife();
        cout << endl;
    }
    cout << endl;
}

/// \brief Sets the limits for the single bateman functions.
void ElementFit::SetSingleBatemanParameterLimits()
{
    //sets limits for N0 of the single functions
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k < i+1; k++)
        {
            //Setting limits for N0
            if(paraVals[k]->GetFixBatemanN0()){
                singleBatemanFitFunctions[i]->FixParameter((k*2), paraVals[k]->GetBatemanN0());
            }else{
                singleBatemanFitFunctions[i]->SetParLimits((k*2), 0., doubleNumEvents*10000);
            }
            //Setting limits for decay constant
            if(paraVals[k]->GetFixDecayConst()){
                singleBatemanFitFunctions[i]->FixParameter((k*2)+1, paraVals[k]->GetDecayConst());
            }else{
                singleBatemanFitFunctions[i]->SetParLimits((k*2)+1, paraVals[k]->GetLowerRangeDecayConst(), paraVals[k]->GetUpperRangeDecayConst());
            }
        }
    }
}

/// \brief Sets the limits for the single integral functions.
void ElementFit::SetSingleIntegralParameterLimits()
{
    //sets limits for N0 of the single functions
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k < i+1; k++)
        {
            //Setting limits for N0
            if(paraVals[k]->GetFixIntegralN0()){
                singleIntegralFitFunctions[i]->FixParameter((k*2), paraVals[k]->GetIntegralN0());
            }else{
                singleIntegralFitFunctions[i]->SetParLimits((k*2), 0., doubleNumEvents*10000);
            }
            //Setting limits for decay constant
            if(paraVals[k]->GetFixDecayConst()){
                singleIntegralFitFunctions[i]->FixParameter((k*2)+1, paraVals[k]->GetDecayConst());
            }else{
                singleIntegralFitFunctions[i]->SetParLimits((k*2)+1, paraVals[k]->GetLowerRangeDecayConst(), (paraVals[k]->GetUpperRangeDecayConst()));
            }
        }
    }
}

/// \brief Sets the limits for the total bateman function.
void ElementFit::SetTotalBatemanParameterLimits()
{
    //DisplayParameterLimits();
    for(int i = 0; i < numElements; i++)
    {
        //1 = case for simulated data, range done with scalar multiples instead of user input range.
        if(inputHistoExecutionType == 1){
            //Setting limits for N0
            if(paraVals[i]->GetFixBatemanN0()){
                batemanFunction->FixParameter((i*2), paraVals[i]->GetBatemanN0());
            }else{
                batemanFunction->SetParLimits((i*2), 0., doubleNumEvents*10000);
            }
            //Setting limits for decay constant
            if(paraVals[i]->GetFixDecayConst()){
                batemanFunction->FixParameter((i*2)+1, paraVals[i]->GetDecayConst());
            }else{
                batemanFunction->SetParLimits((i*2)+1, (paraVals[i]->GetDecayConst() * .01), (paraVals[i]->GetDecayConst() * 100.));
            }
            
        //2 = setting from fit parameters of fitted input histogram, use user defined limits.
        }else if(inputHistoExecutionType == 2 || inputHistoExecutionType == 3)
        {
            //Setting limits for N0
            if(paraVals[i]->GetFixBatemanN0()){
                batemanFunction->FixParameter((i*2), paraVals[i]->GetBatemanN0());
            }else{
                batemanFunction->SetParLimits((i*2), paraVals[i]->GetLowerRangeBatemanN0(), paraVals[i]->GetUpperRangeBatemanN0());
            }
            //Setting limits for decay constant
            if(paraVals[i]->GetFixDecayConst()){
                batemanFunction->FixParameter((i*2)+1, paraVals[i]->GetDecayConst());
            }else{
                batemanFunction->SetParLimits((i*2)+1, (paraVals[i]->GetLowerRangeDecayConst10Ns()), (paraVals[i]->GetUpperRangeDecayConst10Ns()));
            }
        }
    }
}

/// \brief Sets the limits for the total integral function.
void ElementFit::SetTotalIntegralParameterLimits()
{
    for(int i = 0; i < numElements; i++)
    {
        //1 = case for simulated data, range done with scalar multiples instead of user input range.
        if(inputHistoExecutionType == 1){
            //Setting limits for N0
            if(paraVals[i]->GetFixIntegralN0()){
                integralFunction->FixParameter((i*2), paraVals[i]->GetIntegralN0());
            }else{
                integralFunction->SetParLimits((i*2), 0., doubleNumEvents*10000);
            }
            //Setting limits for decay constant
            if(paraVals[i]->GetFixDecayConst()){
                integralFunction->FixParameter((i*2)+1, paraVals[i]->GetDecayConst());
            }else{
                integralFunction->SetParLimits((i*2)+1, (paraVals[i]->GetDecayConst() * .01), (paraVals[i]->GetDecayConst() * 100.));
            }
            
        //2 = setting from fit parameters of fitted input histogram, use user defined limits.
        }else if(inputHistoExecutionType == 2 || inputHistoExecutionType == 3)
        {
            //Setting limits for N0
            if(paraVals[i]->GetFixIntegralN0()){
                integralFunction->FixParameter((i*2), paraVals[i]->GetIntegralN0());
            }else{
                integralFunction->SetParLimits((i*2), paraVals[i]->GetLowerRangeIntegralN0(), paraVals[i]->GetUpperRangeIntegralN0());
            }
            //Setting limits for decay constant
            if(paraVals[i]->GetFixDecayConst()){
                integralFunction->FixParameter((i*2)+1, paraVals[i]->GetDecayConst());
            }else{
                integralFunction->SetParLimits((i*2)+1, (paraVals[i]->GetLowerRangeDecayConst10Ns()), (paraVals[i]->GetUpperRangeDecayConst10Ns()));
            }
        }
    }
}

/// \brief Displays parameter limits set
void ElementFit::DisplayParameterLimits()
{
    cout << "PARAMETER FIT RANGE" << endl;
    //non input histo
    if(inputHistoExecutionType == 1)
    {
        for(int i = 0; i < numElements; i++)
        {
            cout << elementNames[i] << ":" << endl;
            cout << "\tInitial Lower Range: 0" << endl;
            cout << "\tInitial Upper Range: " << doubleNumEvents*10000 << endl;
            cout << "\tHalf Life Lower Range: " << paraVals[i]->GetHalfLife() * .01 << endl;
            cout << "\tHalf Life Upper Range: " << paraVals[i]->GetHalfLife() * 100. << endl;
            cout << endl;
        }
        cout << endl;
    //input histo
    }else if(inputHistoExecutionType == 2 || inputHistoExecutionType == 3)
    {
        for(int i = 0; i < numElements; i++)
        {
            cout << elementNames[i] << ":" << endl;
            cout << "\tBateman initial Lower Range: " << paraVals[i]->GetUpperRangeBatemanN0() << endl;
            cout << "\tBateman initial Upper Range: " << paraVals[i]->GetUpperRangeBatemanN0() << endl;
            cout << "\tIntegral initial Lower Range: " << paraVals[i]->GetUpperRangeIntegralN0() << endl;
            cout << "\tIntegral initial Upper Range: " << paraVals[i]->GetUpperRangeIntegralN0() << endl;
            cout << "\tHalf Life Lower Range: " << paraVals[i]->GetLowerRangeHalfLife() << endl;
            cout << "\tHalf Life Upper Range: " << paraVals[i]->GetUpperRangeHalfLife() << endl;
            cout << endl;
        }
        cout << endl;
    }
}

#endif