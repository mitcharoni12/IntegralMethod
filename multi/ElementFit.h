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
#include "TGraph.h"

#include "SingleElementFitValues.h"
#include "ChainFitValues.h"
#include "SingleCycleHistoHolder.h"
#include "CycleHistoHolder.h"
#include "SingleCycleCanvasHolder.h"
#include "CycleCanvasHolder.h"
#include "SingleCycleGraphHolder.h"
#include "CycleGraphHolder.h"
#include "ParameterValue.h"
#include "FitFunction.h"
#include "FitOption.h"

using namespace std;

typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);

/// This class can be considered the cornerstore for the entire program. This class is the class that generated the data then fits it.
/// After the generation and fitting, the fitted values are passed onto the Run class if multiple Runs/Cycles are chosen.
class ElementFit{
    public:
        ElementFit(Double_t (*batemanFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Double_t (**batemanFitFunctions)(Double_t*, Double_t*), Double_t (**integralFitFunctions)(Double_t*, Double_t*),
                   ParameterValue** paraVals, FitOption* fitOptions);
        ~ElementFit();
        //getter function
        Double_t getElementParameters(int i){return paraVals[i]->getValueDecayConst();}                 ///< Returns decay constant values used to set starting point for fit functions
        FitOption* getFitOptions(){return fitOptions;}                                                  ///< Returns options the user chose to fit with
        ChainFitValues* getBatemanFitParameters(){return totalBatemanFitParameters;}                    ///< Returns fitted values for the fit of the total Bateman histogram
        ChainFitValues* getIntegralFitParameters(){return totalIntegralFitParameters;}                  ///< Returns fitted values for the fit of the total Integral histogram
        SingleElementFitValues* getSingleBatemanFitParameters(){return singleBatemanFitParameters;}     ///< Returns fitted values for the fit of all the single Bateman histograms
        SingleElementFitValues* getSingleIntegralFitParameters(){return singleIntegralFitParameters;}   ///< Returns fitted values for the fit of all the single Integral histograms
        //setter function
        void setNumRuns(Int_t numRuns){this->numRuns = numRuns;}
        void setNumCycles(Int_t numCycles){this->numCycles = numCycles;}
        //public functions
        void displayIntegralGraph(TCanvas* can);
        void displayBatemanHisto(TCanvas* can);
        void displaySingleHistos(TCanvas** can);
        void displayParameters();
        void DrawIndividualHistos(CycleCanvasHolder* batemanTotalCanvases, CycleCanvasHolder* integralTotalCanvases, SingleCycleCanvasHolder* singleBatemanCanvases, SingleCycleCanvasHolder* singleIntegralCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex);
        void fitHistos(Int_t cycleIndex, Int_t runIndex);
        void genAndFillHistos();
    private:
        //private variables
        FitOption* fitOptions;                                                      ///< Contains fit options for program
        string* elementNames;                                                       ///< Contains element names for every element in the decay chain
        Int_t numEvents, numElements, numParameters, numRuns, numCycles;
        bool multiSourceChoice;                                                     ///< True = generate multiple sets of histograms for the different cyles. False = generate a single set of histograms for the diiferent cyles.
        bool rebinChoice;                                                           ///< True = Have program execute with changing bin number between cylces and keeping time fit constant. False = Dont do rebin
        Int_t* binNumArr;                                                           ///< Contains an array of number of bins the histogram have between cycles.
        TF1* integralFunction, *batemanFunction;
        TF1** singleFitFunctions;
        CycleHistoHolder* batemanHisto;                                             ///< Contains every total Bateman histogram for the entire program.
        CycleHistoHolder* integralHisto;                                            ///< Contains every total integral histogram for the entire program.
        SingleCycleHistoHolder* singleBatemanHisto;                                 ///< Contains every single Bateman histogram for the entire program.
        SingleCycleHistoHolder* singleIntegralHisto;                                ///< Contains every single integral histogram for the entire program.
        CycleGraphHolder* integralGraph;                                            ///< Contains every total integral graph for the entire program.
        SingleCycleGraphHolder* singleIntegralGraph;                                ///< Contains every single integral graph for the entire program.
        //passed fit functions used to make the TF1 for fitting
        decayFunction* batemanFitFunctions, *integralFitFunctions;
        decayFunction passedBatemanFunction, passedIntegralFunction;
        TRandom3 rand;                                                              ///< Used for random number generation.
        Double_t timeRunStart;                                                      ///< Inital value for the start of fit time, can change between cycle.
        Double_t timeRunEnd;                                                        ///< Inital value for the end of fit time, can change between cycle.
        Double_t leaveOutStartBinNumber;                                            ///< How many bins to leave out of the start of the fit.
        Double_t leaveOutEndBinNumber;                                              ///< How many bins to leave out of the end of the fit.
        Double_t doubleNumEvents;                                                   ///< Number of events except its of data type Double_t.
        Double_t* randArr;                                                          ///< Array used for storing random values for event generation.
        Double_t* binWidth;                                                         ///< Holds bin widths between cycles.
        Double_t* timeEndArr;                                                       ///< Holds the time of fit end between cycles.
        Double_t* timeStartArr;                                                     ///< Holds the time of fit start between cycles.
        ChainFitValues* totalBatemanFitParameters, *totalIntegralFitParameters;
        SingleElementFitValues* singleBatemanFitParameters, *singleIntegralFitParameters;
        ParameterValue** paraVals;                                                  ///< Accepted values for all parameters in all the fit functions.
        Int_t globalSeedChanger = 0;                                                ///< Used for assigning in changing the seed.
        //helper functions
        void changeSeed();
        void createIntegralGraph();
        void createHistoHolders();
        void createSingleFitFunctions(Int_t timeEnd);
        void createTotalFitFunctions(Int_t timeEnd);
        void fitIntegralHisto(Int_t cycleIndex, Int_t runIndex, Double_t startFit, Double_t endFit);
        void fitBatemanHisto(Int_t cycleIndex, Int_t runIndex, Double_t startFit, Double_t endFit);
        void fitSingleHistos(Int_t cycleIndex, Int_t runIndex, Double_t startFit, Double_t endFits);
        void DisplayParameterLimits();
        void DisplayTotalFunctionParameters();
        void genBatemanHistograms();
        void genIntegralHisto();
        void genIntegralSingleHistos();
        void genIntegralHistoSimulated();
        void setFunctionParametersTotal();
        void setFunctionParamersSingle();
        void setParaLimits();
        //TCanvas* test = new TCanvas("test", "test", 500, 500);
        clock_t Tclock;
};  

/// \brief The constructor
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
    this->timeRunEnd = fitOptions->GetTimeRunEnd();
    this->paraVals = paraVals;
    this->multiSourceChoice = fitOptions->GetMultiSourceChoice();
    this->rebinChoice = fitOptions->GetRebinChoice();
    this->leaveOutStartBinNumber = fitOptions->GetStartBinLeaveOut();
    this->leaveOutEndBinNumber = fitOptions->GetEndBinLeaveOut();
    doubleNumEvents = (Double_t) numEvents;

    //creating fit parameter holders
    singleBatemanFitParameters = new SingleElementFitValues(numElements);
    singleIntegralFitParameters = new SingleElementFitValues(numElements);
    totalBatemanFitParameters = new ChainFitValues(numElements);
    totalIntegralFitParameters = new ChainFitValues(numElements);

    singleFitFunctions = new TF1* [numElements*2];
    //setting values for randomization
    randArr = new Double_t [numElements];
    Tclock = clock();
    //creates required arrays for number of bins and time fit
    fitOptions->CreateRequiredDataSets();
    binWidth = fitOptions->GetBinWidthArr();
    timeEndArr = fitOptions->GetTimeFitEndArr();
    timeStartArr = fitOptions->GetTimeFitStartArr();
    binNumArr = fitOptions->GetBinNumArr();
    //generates all the histo objects and 
    genAndFillHistos();
    //removes restrictions on fitting function calls and itterations
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000000);
}

ElementFit::~ElementFit()
{
    delete singleBatemanFitParameters;
    delete singleIntegralFitParameters;
    delete totalBatemanFitParameters;
    delete totalIntegralFitParameters;
    for(int i = 0; i < numElements*2; i++)
    {
        delete singleFitFunctions[i];
    }
    delete [] singleFitFunctions;
    delete [] randArr;
    delete batemanFunction;
    delete integralFunction;
    delete batemanHisto;
    delete integralHisto;
    delete singleBatemanHisto;
    delete singleIntegralHisto;
    delete integralGraph;
    delete singleIntegralGraph;
}

/// \brief Creates required histograms for the ENTIRE program.
///
/// Creates all the histogram objects that will be required for the program. When doing things like changing binning and time run the histograms are required to be different.
void ElementFit::createHistoHolders()
{
    string histoName;
    histoName = "Total Bateman Histo";
    batemanHisto = new CycleHistoHolder(numCycles, numRuns, histoName, binNumArr, timeEndArr);
    histoName = "Total Integral Histo";
    integralHisto = new CycleHistoHolder(numCycles, numRuns, histoName, binNumArr, timeEndArr);
    histoName = "Single Bateman Histo";
    singleBatemanHisto = new SingleCycleHistoHolder(numCycles, numElements, numRuns, histoName, binNumArr, timeEndArr, elementNames);
    histoName = "Single Integral Histo";
    singleIntegralHisto = new SingleCycleHistoHolder(numCycles, numElements, numRuns, histoName, binNumArr, timeEndArr, elementNames);
}

/// \brief creates the integral graphs from integral histogram
///
/// creates integral graphs from integral histogram. Does this for all integral histograms once function is called.
void ElementFit::createIntegralGraph()
{
    TH1D* tempHisto;
    TGraph* tempGraph;
    Int_t* pointsArr = new Int_t [numCycles];
    Double_t tempBinWidth, tempTimeValue, tempBinValue;
    for(int i = 0; i < numCycles; i++)
    {
        pointsArr[i] = binNumArr[i] + 1;
    }

    //creates the storage objects for in the integral graphs
    integralGraph = new CycleGraphHolder(numCycles, numRuns, "Total Integral Graph", pointsArr);
    singleIntegralGraph = new SingleCycleGraphHolder(numCycles, numElements, numRuns, "Single Integral Graph", pointsArr, elementNames);

    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        //bin width for cycle
        tempBinWidth = binWidth[cycleIndex];

        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {  
            //storing graph object that will be created from the histogram
            tempHisto = integralHisto->GetAHisto(cycleIndex, runIndex);
            tempGraph = integralGraph->GetAGraph(cycleIndex, runIndex);
            tempGraph->SetPoint(1, 0.0, 0.0);
            
            //creates graph for total integral graph
            for(int i = 1; i < binNumArr[cycleIndex] + 1; i++)
            {
                tempTimeValue = i * tempBinWidth;
                tempBinValue = tempHisto->GetBinContent(i);
                tempGraph->SetPoint(i+1, tempTimeValue, tempBinValue);
            }
            //creates the graphs for the single integral graphs
            for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
            {
                tempHisto = singleIntegralHisto->GetAHisto(cycleIndex, runIndex, elementIndex);
                tempGraph = singleIntegralGraph->GetAGraph(cycleIndex, runIndex, elementIndex);
                tempGraph->SetPoint(1, 0.0, 0.0);

                for(int i = 1; i < binNumArr[cycleIndex] + 1; i++)
                {
                    tempTimeValue = i * tempBinWidth;
                    tempBinValue = tempHisto->GetBinContent(i);
                    tempGraph->SetPoint(i+1, tempTimeValue, tempBinValue);
                }
            }
        }
    }

    delete [] pointsArr;
}

/// \brief changes seed for random number generator
void ElementFit::changeSeed()
{
    Double_t seederSeed;
    Double_t randNumber;

    randNumber = rand.Uniform();
    seederSeed = Tclock + globalSeedChanger + randNumber;

    rand.SetSeed(Tclock + globalSeedChanger + seederSeed);
    globalSeedChanger++;
}

/// \brief creates the TF1 objects for the single Bateman and integral fit functions 
/// creates the TF1 objects for the single Bateman and integral fit functions. Creates these from the passed fit function.
/// Must create these because you have to fit with a TF1 function.
void ElementFit::createSingleFitFunctions(Int_t timeEnd)
{
    for(int i = 0; i < numElements; i++)
    {
        singleFitFunctions[(i*2)] = new TF1((elementNames[i] + "RegSingFunc").c_str(), batemanFitFunctions[i], 0., timeEnd, (i+1)*2);
        singleFitFunctions[(i*2)+1] = new TF1((elementNames[i] + "InteSingFunc").c_str(), integralFitFunctions[i], 0., timeEnd, (i+1)*2);
    }
}

/// \brief creates the TF1 objects for the total Bateman and integral fit functions 
/// creates the TF1 objects for the total Bateman and integral fit functions. Creates these from the passed fit function.
/// Must create these because you have to fit with a TF1 function.
void ElementFit::createTotalFitFunctions(Int_t timeEnd)
{
    batemanFunction = new TF1("TotalBatemanFunction", passedBatemanFunction, 0., timeEnd, numParameters);
    integralFunction = new TF1("TotalIntegralFunction", passedIntegralFunction, 0., timeEnd, numParameters);
}

/// \brief Displays total integral histogram for the single fit option for the program
///
/// Displays the total integral graph for single run program execution type
void ElementFit::displayIntegralGraph(TCanvas* can)
{
    can->cd();
    integralGraph->GetAGraph(0, 0)->Draw();
}

/// \brief Displays fitted parameters for the single fit option for the program
///
/// Displays fitted parameters for the single run program execution type
void ElementFit::displayParameters()
{  
    cout << endl << "BATEMAN FIT PARAMETERS/ERRORS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ": " << "\tHalfLife: " << totalBatemanFitParameters->GetAnHalfLife(i) << "s" << endl << 
        "\tError: " << totalBatemanFitParameters->GetAnHalfLifeError(i) << "s" << endl;
        cout << "\tN0: " << totalBatemanFitParameters->GetAnN0(i) << "" << endl << 
        "\tError: " << totalBatemanFitParameters->GetAnN0Error(i) << "" << endl << endl;
    }
    cout << "INTEGRAL FIT PARAMETERS/ERRORS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ": " << "\tHalfLife: " << totalIntegralFitParameters->GetAnHalfLife(i) << "s" << endl << 
        "\tError: " << totalIntegralFitParameters->GetAnHalfLifeError(i) << "s" << endl;
        cout << "\tN0: " << totalIntegralFitParameters->GetAnN0(i) << "" << endl << 
        "\tError: " << totalIntegralFitParameters->GetAnN0Error(i) << "" << endl << endl;
    }
    cout << "BATEMAN SINGLE FIT PARAMETERS/ERRORS" << endl << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << " Single Fit Values" << endl;
        for(int k = 0; k < i+1; k++)
        {
            cout << elementNames[k] << ": \tHalf Life: " << singleBatemanFitParameters->GetAnHalfLife(i, k) << "s" << endl << 
            "\tHalf Life Error: " << singleBatemanFitParameters->GetAnHalfLifeError(i,k) << "s" << endl;
            cout << "\tN0: " << singleBatemanFitParameters->GetAnN0(i, k) << "" << endl << 
            "\tN0 Error: " << singleBatemanFitParameters->GetAnN0Error(i,k) << "" << endl;
            cout << endl;
        }
    }
    cout << endl << "INTEGRAL SINGLE FIT PARAMETERS/ERRORS" << endl << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << " Single Fit Values" << endl;
        for(int k = 0; k < i+1; k++)
        {
            cout << elementNames[k] << ": \tHalf Life: " << singleIntegralFitParameters->GetAnHalfLife(i, k) << "s" << endl << 
            "\tHalf Life Error: " << singleIntegralFitParameters->GetAnHalfLifeError(i,k) << "s" << endl;
            cout << "\tN0: " << singleIntegralFitParameters->GetAnN0(i, k) << "" << endl << 
            "\tN0 Error: " << singleIntegralFitParameters->GetAnN0Error(i,k) << "" << endl;
            cout << endl;
        }
    }
}

/// \brief Displays total Bateman histogram for the single fit option for the program
///
/// Displays the total Bateman histogram for single run program execution type
void ElementFit::displayBatemanHisto(TCanvas* can)
{
    can->cd();
    batemanHisto->GetAHisto(0, 0)->Draw();
}

/// \brief Displays single histograms for the single fit option for the program
///
/// Displays the single Bateman and integral histogram/graphs for single run program execution type
void ElementFit::displaySingleHistos(TCanvas** can)
{
    for(int i = 0; i < numElements; i++)
    {
        can[(i*2)]->cd();
        singleBatemanHisto->GetAHisto(0, 0, i)->Draw();

        can[(i*2)+1]->cd();
        singleIntegralGraph->GetAGraph(0, 0, i)->Draw();
    }
}

/// \brief Displays the individual histogram for any execution type of the program
///
/// Draws the invidivual histograms for the program within the range set by the user
void ElementFit::DrawIndividualHistos(CycleCanvasHolder* batemanTotalCanvases, CycleCanvasHolder* integralTotalCanvases, SingleCycleCanvasHolder* singleBatemanCanvases, SingleCycleCanvasHolder* singleIntegralCanvases,
                                      Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex)
{
    upperRunIndex = upperRunIndex + 1;
    upperCycleIndex = upperCycleIndex + 1;

    for(int cycleIndex = lowerCycleIndex; cycleIndex < upperCycleIndex; cycleIndex++)
    {
        for(int runIndex = lowerRunIndex; runIndex < upperRunIndex; runIndex++)
        {
            batemanTotalCanvases->GetACanvas(cycleIndex, runIndex)->cd();
            batemanHisto->GetAHisto(cycleIndex, runIndex)->Draw();

            integralTotalCanvases->GetACanvas(cycleIndex, runIndex)->cd();
            integralGraph->GetAGraph(cycleIndex, runIndex)->Draw();

            for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
            {
                singleBatemanCanvases->GetACanvas(cycleIndex, runIndex, elementIndex)->cd();
                singleBatemanHisto->GetAHisto(cycleIndex, runIndex, elementIndex)->Draw();

                singleIntegralCanvases->GetACanvas(cycleIndex, runIndex, elementIndex)->cd();
                singleIntegralGraph->GetAGraph(cycleIndex, runIndex, elementIndex)->Draw();
            }
        }
    }
}

/// \brief Calls fit methods to fit histograms at specific run and cycle index
///
/// Fits the single/total Bateman and integral histograms/graphs for a specified cycle and run index. The logic for leaving out certain bins in also included here.
/// The fit function parameters are also reset before the fit here
void ElementFit::fitHistos(Int_t cycleIndex, Int_t runIndex)
{
    timeRunStart = timeStartArr[cycleIndex];
    timeRunEnd = timeEndArr[cycleIndex];
    createSingleFitFunctions(timeEndArr[cycleIndex]);
    createTotalFitFunctions(timeEndArr[cycleIndex]);
    setFunctionParametersTotal();
    setFunctionParamersSingle();
    setParaLimits();

    Double_t startFitOffset;
    Double_t startFit;
    startFitOffset = leaveOutStartBinNumber * binWidth[cycleIndex];
    startFit = startFitOffset + timeRunStart;
    Double_t endFitOffset;
    Double_t endFit;
    endFitOffset = leaveOutEndBinNumber * binWidth[cycleIndex];
    endFit = timeRunEnd - endFitOffset;

    cout << "CYCLE: " << cycleIndex << " START FIT: " << startFit << " END FIT: " << endFit << endl;

    fitBatemanHisto(cycleIndex, runIndex, startFit, endFit);
    fitIntegralHisto(cycleIndex, runIndex, startFit, endFit);
    fitSingleHistos(cycleIndex, runIndex, startFit, endFit);
}

/// \brief Fits total integral histogram at specified cycle and run index
///
/// Fits the total integral histogram at a specified cycle and run index.
/// The fit parameters are extracted afer the fit.
void ElementFit::fitIntegralHisto(Int_t cycleIndex, Int_t runIndex, Double_t startFit, Double_t endFit)
{
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;
    TH1D* tempHisto;
    TGraph* tempGraph;

    cout << "FITTING TOTAL INTEGRAL CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;

    tempGraph = integralGraph->GetAGraph(cycleIndex, runIndex);
    tempGraph->Fit(integralFunction, "", "", startFit, endFit);

    for(int i = 0; i < numElements; i++)
    {   
        valueN0 = integralFunction->GetParameter((i*2));
        errorN0 = integralFunction->GetParError((i*2));
        valueDecayConst = integralFunction->GetParameter((i*2)+1);
        errorDecayConst = integralFunction->GetParError((i*2)+1);
        valueHalfLife = log(2)/(valueDecayConst);
        errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));
        //used to get either the N0 value or half life value
        totalIntegralFitParameters->SetAnN0(i, valueN0);
        totalIntegralFitParameters->SetAnN0Error(i, errorN0);
        totalIntegralFitParameters->SetAnHalfLife(i, valueHalfLife);
        totalIntegralFitParameters->SetAnHalfLifeError(i, errorHalfLife);
    }
}

/// \brief Fits total Bateman histogram at specified cycle and run index
///
/// Fits total Bateman histogram at specified cycle and run index.
/// The fit parameters are extracted afer the fit.
void ElementFit::fitBatemanHisto(Int_t cycleIndex, Int_t runIndex, Double_t startFit, Double_t endFit)
{
    //dynamic array
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;
    TH1D* tempHisto;

    cout << "FITTING TOTAL BATEMAN CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;

    tempHisto = batemanHisto->GetAHisto(cycleIndex, runIndex);
    tempHisto->Fit(this->batemanFunction, "L", "", timeRunStart, endFit);
    for(int i = 0; i < numElements; i++)
    {   
        valueN0 = batemanFunction->GetParameter((i*2));
        errorN0 = batemanFunction->GetParError((i*2));
        valueDecayConst = batemanFunction->GetParameter((i*2)+1);
        errorDecayConst = batemanFunction->GetParError((i*2)+1);
        valueHalfLife = log(2)/(valueDecayConst);
        errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));
        //used to get either the N0 value or half life value
        totalBatemanFitParameters->SetAnN0(i, valueN0);
        totalBatemanFitParameters->SetAnN0Error(i, errorN0);
        totalBatemanFitParameters->SetAnHalfLife(i, valueHalfLife);
        totalBatemanFitParameters->SetAnHalfLifeError(i, errorHalfLife);
    }
}

/// \brief Fits single Histograms at specified cycle and run index
///
/// Fits single Histograms at specified cycle and run index.
/// The fit parameters are extracted afer the fit.
void ElementFit::fitSingleHistos(Int_t cycleIndex, Int_t runIndex, Double_t startFit, Double_t endFit)
{
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;
    TH1D* tempSingleBatemanHisto, *tempSingleIntegralHisto;
    TGraph* tempSingleGraph;

    //loop for fitting and storing value for elements
    for(int i = 0; i < numElements; i++)
    {
        tempSingleBatemanHisto = singleBatemanHisto->GetAHisto(cycleIndex, runIndex, i);
        tempSingleGraph = singleIntegralGraph->GetAGraph(cycleIndex, runIndex, i);

        cout << "FITTING SINGLE BATEMAN " << elementNames[i] << " CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;
        tempSingleBatemanHisto->Fit(singleFitFunctions[(i*2)], "L", "", startFit, endFit);
        cout << "FITTING SINGLE INTEGRAL " << elementNames[i] << " CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;
        tempSingleGraph->Fit(singleFitFunctions[(i*2)+1], "", "", startFit, endFit);

        for(int k = 0; k < (i+1); k++)
        {
            //storing values for bateman part of function
            valueN0 = singleFitFunctions[(i*2)]->GetParameter((k*2));
            errorN0 = singleFitFunctions[(i*2)]->GetParError((k*2));
            valueDecayConst = singleFitFunctions[(i*2)]->GetParameter((k*2)+1);
            errorDecayConst = singleFitFunctions[(i*2)]->GetParError((k*2)+1);
            valueHalfLife = log(2)/(valueDecayConst);
            errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));

            singleBatemanFitParameters->SetAnN0(i, k, valueN0);
            singleBatemanFitParameters->SetAnN0Error(i, k, errorN0);
            singleBatemanFitParameters->SetAnHalfLife(i, k, valueHalfLife);
            singleBatemanFitParameters->SetAnHalfLifeError(i, k, errorHalfLife);

            //storing values for integral part of function
            valueN0 = singleFitFunctions[(i*2)+1]->GetParameter((k*2));
            errorN0 = singleFitFunctions[(i*2)+1]->GetParError((k*2));
            valueDecayConst = singleFitFunctions[(i*2)+1]->GetParameter((k*2)+1);
            errorDecayConst = singleFitFunctions[(i*2)+1]->GetParError((k*2)+1);
            valueHalfLife = log(2)/(valueDecayConst);
            errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));

            singleIntegralFitParameters->SetAnN0(i, k, valueN0);
            singleIntegralFitParameters->SetAnN0Error(i, k, errorN0);
            singleIntegralFitParameters->SetAnHalfLife(i, k, valueHalfLife);
            singleIntegralFitParameters->SetAnHalfLifeError(i, k, errorHalfLife);
        }
    }
}

/// \brief Function to call all the histograms in the program to be generated
///
/// Function to call all the histograms in the program to be generated.
/// Just calls to sub functions, real code are in the indivudal sub functions.
void ElementFit::genAndFillHistos()
{
    createHistoHolders();
    genBatemanHistograms();
    genIntegralHistoSimulated();
    createIntegralGraph();
}

/// \brief Generates all total integral histograms for total Bateman histograms
///
/// Generates all total integral histograms for total Bateman histograms.
/// Each individual bin(i) in the integral histogram is the sum of the integral bin i-1 plus the bin i from the regular histogram.
void ElementFit::genIntegralHisto()
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

            tempIntegralHisto->SetBinContent(1, tempBatemanHisto->GetBinContent(1));
            for(int i = 2; i<= tempNumBins; i++)
            {
                tempIntegralHisto->SetBinContent(i, tempBatemanHisto->GetBinContent(i) + tempIntegralHisto->GetBinContent(i-1));
            }
        }
    }
}

/// \brief Generates all integral histograms from the Bateman histograms
///
/// Generates all integral histograms from the Bateman histograms.
/// Just call to sub functions.
void ElementFit::genIntegralHistoSimulated()
{
    genIntegralHisto();
    genIntegralSingleHistos();
}

/// \brief Generates all the single integral histograms from the single Bateman histograms.
///
/// Generates all the single integral histograms from the single Bateman histograms.
/// Each individual bin(i) in the integral histogram is the sum of the integral bin i-1 plus the bin i from the regular histogram.
void ElementFit::genIntegralSingleHistos()
{
    TH1D* tempIntegralHisto, *tempBatemanHisto;
    Int_t tempNumBins;

    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        tempNumBins = binNumArr[cycleIndex];
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            for(int i = 0; i < numElements; i++)
            {   
                tempIntegralHisto = singleIntegralHisto->GetAHisto(cycleIndex, runIndex, i);
                tempBatemanHisto = singleBatemanHisto->GetAHisto(cycleIndex, runIndex, i);

                tempIntegralHisto->SetBinContent(1, tempBatemanHisto->GetBinContent(1));
                for(int k = 2; k <= tempNumBins; k++)
                {
                    tempIntegralHisto->SetBinContent(k, tempIntegralHisto->GetBinContent(k-1) + tempBatemanHisto->GetBinContent(k));
                }
            }
        }
    }
}

/// \brief Generates the total and single Bateman histograms
///
/// Generates the total and single Bateman histograms. Has different methods of generation for different program execution types.
/// Still not sure why generation works. This is the single most important function in the entire program.
void ElementFit::genBatemanHistograms()
{
    TH1D* tempHisto;
    TH1D** singleTempHisto;
    singleTempHisto = new TH1D* [numElements];
    Double_t hold = 0;
    Double_t stack = 0;
    //case for generating the single histogram for the single source histogram choice
    if(!rebinChoice && !multiSourceChoice)
    {
        TH1D* tempHisto;
        TH1D** tempSingleHisto;
        tempSingleHisto = new TH1D* [numElements];

        changeSeed();

        //generates events for the single and total histograms symotaniously
        for(int i = 0; i < numEvents; i++)
        {
            //generating the times the events occured
            for(int j = 0; j < numElements; j++)
            {
                hold = rand.Uniform();
                randArr[j] = (-TMath::Log(hold)) / (paraVals[j]->getValueDecayConst());
            }
            
            //putting events in repsective histograms
            for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
            {
                stack += randArr[elementIndex];
                for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
                {
                    tempHisto = batemanHisto->GetAHisto(cycleIndex, 0);
                    singleTempHisto[elementIndex] = singleBatemanHisto->GetAHisto(cycleIndex, 0, elementIndex);
                    singleTempHisto[elementIndex]->Fill(stack);
                    tempHisto->Fill(stack);
                }
            }
            stack = 0.0f; 
        }

        delete [] tempSingleHisto;
    //case for multiple histogram for the multiple source histogram choice.
    }else if(!rebinChoice && multiSourceChoice)
    {
        for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
        {
            for(int runIndex = 0; runIndex < numRuns; runIndex++)
            {
                //getting histograms to generate events for
                tempHisto = batemanHisto->GetAHisto(cycleIndex, runIndex);
                for(int i = 0; i < numElements; i++)
                {
                    singleTempHisto[i] = singleBatemanHisto->GetAHisto(cycleIndex, runIndex, i);
                }

                changeSeed();
                //generates events for the single and total histograms symotaniously
                for(int i = 0; i < numEvents; i++)
                {
                    //generating the times the events occured
                    for(int j = 0; j < numElements; j++)
                    {
                        hold = rand.Uniform();
                        randArr[j] = (-TMath::Log(hold)) / (paraVals[j]->getValueDecayConst());
                    }

                    //putting events in repsective histograms
                    for(int k = 0; k < numElements; k++)
                    {
                        stack += randArr[k];
                        singleTempHisto[k]->Fill(stack);
                        tempHisto->Fill(stack); 
                    }
                    stack = 0.0f;
                }
            }
        }
    //case for generating data with rebinning
    }else if(rebinChoice)
    {
        //numbers generated same but data fed in differently. We want all the events in every run of each cycle to be identical so we can see the effects of rebinning. 
        //therefore data generated in run index 0 of cycle 0 must be identical to run index 0 of cycle 1 and so on.
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            //change seed for generation between each run
            changeSeed();
            //generates events for the single and total histograms symotaniously
            for(int i = 0; i < numEvents; i++)
            {
                //generating the times the events occured
                for(int j = 0; j < numElements; j++)
                {
                    hold = rand.Uniform();
                    randArr[j] = (-TMath::Log(hold)) / (paraVals[j]->getValueDecayConst());
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
                        singleTempHisto[k] = singleBatemanHisto->GetAHisto(cycleIndex, runIndex, k);
                        singleTempHisto[k]->Fill(stack);
                        tempHisto->Fill(stack);
                    }
                }
                stack = 0.0f;
            }
        }
    }

    delete [] singleTempHisto;
}

/// \brief Sets the parameters for the total fit functions
///
/// Sets the parameters for the total fit functions. Must do this with a TF1 or else it will not do the fit. Just set the initial values to the accepted values.
void ElementFit::setFunctionParametersTotal()
{
    DisplayTotalFunctionParameters();
    for(int i = 0; i < numElements; i++)
    {
        //have to have the if statments to account for the fact that the value could be a range
        if(paraVals[i]->getIsValueInitValue())
        {
            integralFunction->SetParameter((i*2), paraVals[i]->getInitValue());
            batemanFunction->SetParameter((i*2), paraVals[i]->getInitValue());
        }else{
            integralFunction->SetParameter((i*2), paraVals[i]->getRangeAverageInitValue());
            batemanFunction->SetParameter((i*2), paraVals[i]->getRangeAverageInitValue());
        }
        if(paraVals[i]->getIsValueDecayConst())
        {
            integralFunction->SetParameter((i*2)+1, paraVals[i]->getValueDecayConst());
            batemanFunction->SetParameter((i*2)+1, paraVals[i]->getValueDecayConst());
        }else{
            integralFunction->SetParameter((i*2)+1, paraVals[i]->getRangeAverageDecayConst());
            batemanFunction->SetParameter((i*2)+1, paraVals[i]->getRangeAverageDecayConst());
        }
    }
}

/// \brief Displays the parameters used to set the fit functions
void ElementFit::DisplayTotalFunctionParameters()
{
    cout << "INITIAL FIT FUNCTION PARAMETERS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ":\tInit value: " << paraVals[i]->getInitValue();
        cout << "\tHalf Life: " << paraVals[i]->getValueHalfLife();
        cout << endl;
    }
    cout << endl;
}

/// \brief Sets the parameters for the single fit functions
///
/// Sets the parameters for the single fit functions. Must do this with a TF1 or else it will not do the fit. Just set the initial values to the accepted values.
void ElementFit::setFunctionParamersSingle()
{
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k <= i; k++)
        {
            //need if statments to account for the fact that the value could be range
            if(paraVals[k]->getIsValueInitValue())
            {
                singleFitFunctions[(i*2)]->SetParameter((k*2), paraVals[k]->getInitValue());
                singleFitFunctions[(i*2)+1]->SetParameter((k*2), paraVals[k]->getInitValue());
            }else{
                singleFitFunctions[(i*2)]->SetParameter((k*2), paraVals[k]->getRangeAverageInitValue());
                singleFitFunctions[(i*2)+1]->SetParameter((k*2), paraVals[k]->getRangeAverageInitValue());
            }
            if(paraVals[k]->getIsValueDecayConst())
            {
                singleFitFunctions[(i*2)]->SetParameter(((k*2)+1), paraVals[k]->getValueDecayConst());
                singleFitFunctions[(i*2)+1]->SetParameter(((k*2)+1), paraVals[k]->getValueDecayConst());
            }else{
                singleFitFunctions[(i*2)]->SetParameter(((k*2)+1), paraVals[k]->getRangeAverageDecayConst());
                singleFitFunctions[(i*2)+1]->SetParameter(((k*2)+1), paraVals[k]->getRangeAverageDecayConst());
            }
        }
    }
}

/// \brief Sets the limits in which the parameters can fit in the fit function.
///
/// Sets the limits in which the parameters can fit in the fit function. If they are too large the program will not fit correctly but too large then its like cheating.
/// I just use scalar multiples of the accepted values in the range.
void ElementFit::setParaLimits()
{
    DisplayParameterLimits();
    //sets limits for N0 of the total function
    for(int i = 0; i < numElements; i++)
    {
        if(paraVals[i]->getIsValueInitValue())
        {
            batemanFunction->SetParLimits((i*2), 0., doubleNumEvents*10000);
            integralFunction->SetParLimits((i*2), 0., doubleNumEvents*10000);
        }else{
            batemanFunction->SetParLimits((i*2), paraVals[i]->getLowerRangeInitValue(), paraVals[i]->getUpperRangeInitValue());
            integralFunction->SetParLimits((i*2), paraVals[i]->getLowerRangeInitValue(), paraVals[i]->getUpperRangeInitValue());
        }
    }
    //sets limits for lambda value of total function
    for(int i = 0; i < numElements; i++)
    {
        if(paraVals[i]->getIsValueDecayConst())
        {
            batemanFunction->SetParLimits((i*2)+1, (paraVals[i]->getValueDecayConst() * .01), (paraVals[i]->getValueDecayConst() * 100.));
            integralFunction->SetParLimits((i*2)+1, (paraVals[i]->getValueDecayConst() * .01), (paraVals[i]->getValueDecayConst() * 100.));
        }else{
            batemanFunction->SetParLimits((i*2)+1, (paraVals[i]->getLowerRangeDecayConst()), (paraVals[i]->getUpperRangeDecayConst()));
            integralFunction->SetParLimits((i*2)+1, (paraVals[i]->getLowerRangeDecayConst()), (paraVals[i]->getUpperRangeDecayConst()));
        }
    }

    //sets limits for N0 of the single functions
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k < i+1; k++)
        {
            if(paraVals[i]->getIsValueInitValue())
            {
                (singleFitFunctions[(i*2)])->SetParLimits((k*2), 0., doubleNumEvents*10000);
                (singleFitFunctions[(i*2)+1])->SetParLimits((k*2), 0., doubleNumEvents*10000);
            }else{
                (singleFitFunctions[(i*2)])->SetParLimits((k*2), paraVals[i]->getLowerRangeInitValue(), paraVals[i]->getUpperRangeInitValue());
                (singleFitFunctions[(i*2)+1])->SetParLimits((k*2), paraVals[i]->getLowerRangeInitValue(), paraVals[i]->getUpperRangeInitValue());
            }
        }
    }
    //sets limits for lambda values of the single functions
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k <= i; k++)
        {
            if(paraVals[i]->getIsValueDecayConst())
            {
                (singleFitFunctions[(i*2)])->SetParLimits((k*2)+1, (paraVals[k]->getValueDecayConst() * .01), (paraVals[k]->getValueDecayConst() * 100.));
                (singleFitFunctions[(i*2)+1])->SetParLimits((k*2)+1, (paraVals[k]->getValueDecayConst() * .01), (paraVals[k]->getValueDecayConst() * 100.));
            }else{
                (singleFitFunctions[(i*2)])->SetParLimits((k*2)+1, (paraVals[i]->getLowerRangeDecayConst()), (paraVals[i]->getUpperRangeDecayConst()));
                (singleFitFunctions[(i*2)+1])->SetParLimits((k*2)+1, (paraVals[i]->getLowerRangeDecayConst()), (paraVals[i]->getUpperRangeDecayConst()));
            }
        }
    }
}

/// \brief Displays parameter limits set
void ElementFit::DisplayParameterLimits()
{
    cout << "PARAMETER FIT RANGE" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ":" << endl;
        cout << "\tInitial Lower Range: 0" << endl;
        cout << "\tInitial Lower Range: " << doubleNumEvents*10000 << endl;
        cout << "\tHalf Life Lower Range: " << paraVals[i]->getValueHalfLife() * .01 << endl;
        cout << "\tHalf Life Upper Range: " << paraVals[i]->getValueHalfLife() * 100. << endl;
        cout << endl;
    }
}

#endif