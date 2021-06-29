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
#include "SingleCycleHistoHolder.h"
#include "ParameterValue.h"
#include "ChainFitValues.h"
#include "FitValues.h"
#include "SingleElementFitValues.h"
#include "CycleHistoHolder.h"
#include "CycleCanvasHolder.h"
#include "SingleCycleCanvasHolder.h"

using namespace std;

typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);

class ElementFit{
    public:
        ElementFit(Int_t events, Int_t numRuns, Int_t numCycles, Double_t (*regularFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Double_t (**fitFunctions)(Double_t*, Double_t*),
                   Int_t numElements, Double_t timeRunEnd, Int_t numBins, string* elementNames,  ParameterValue** paraVals, Int_t singleHistoChoice, Int_t rebinChoice, Int_t rebinDifference);
        ~ElementFit();
        //getter function
        Double_t getElementParameters(int i){return paraVals[i]->getValueDecayConst();}
        Int_t getNumEvents(){return events;}
        Int_t getNumElements(){return numElements;}
        Int_t* getBinArr(){return binSizeArr;}
        Int_t getTimeRunEnd(){return timeRunEnd;}
        ChainFitValues* getRegularFitParameters(){return totalRegularFitParameters;}
        ChainFitValues* getIntegralFitParameters(){return totalIntegralFitParameters;}
        SingleElementFitValues* getSingleRegularFitParameters(){return singleRegularFitParameters;}
        SingleElementFitValues* getSingleIntegralFitParameters(){return singleIntegralFitParameters;}
        //setter function
        //void setNumBins(Int_t numBins);
        void setTimeRunStart(Double_t timeRunStart){this->timeRunStart = timeRunStart;}
        void setTimeRunEnd(Double_t timeRunEnd){this->timeRunEnd = timeRunEnd;}
        void setNumRuns(Int_t numRuns){this->numRuns = numRuns;}
        void setNumCycles(Int_t numCycles){this->numCycles = numCycles;}
        //public functions
        void createHistoHolders();
        void displayIntegralHisto(TCanvas* can);
        void displayRegularHisto(TCanvas* can);
        void displaySingleHistos(TCanvas** can);
        void displayParameters();
        void DrawIndividualHistos(CycleCanvasHolder* regularTotalCanvases, CycleCanvasHolder* integralTotalCanvases, SingleCycleCanvasHolder* singleRegularCanvases, SingleCycleCanvasHolder* singleIntegralCanvases, Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex);
        void fitDataGenOnce(Int_t cycleIndex, Int_t runIndex);
        void fitHistos(Int_t cycleIndex, Int_t runIndex);
        void fitIntegralHisto(Int_t cycleIndex, Int_t runIndex, Double_t leaveOutStartBinNumber, Double_t leaveOutEndBinNumber);
        void fitRegularHisto(Int_t cycleIndex, Int_t runIndex, Double_t leaveOutStartBinNumber, Double_t leaveOutEndBinNumber);
        void fitSingleHistos(Int_t cycleIndex, Int_t runIndex, Double_t leaveOutStartBinNumber, Double_t leaveOutEndBinNumber);
        void genIntegralHisto();
        void genIntegralSingleHistos();
        void genRandomAlternate();
        void genAndFillHistos();
    private:
        //private variables
        string* elementNames;
        Int_t events, numElements, numBins, numParameters, numRuns, numCycles, singleHistoChoice, rebinChoice, rebinDifference;
        Int_t* binSizeArr;
        TF1* integralFunction, *regularFunction;
        TF1** singleFitFunctions;
        CycleHistoHolder* regularHisto, *integralHisto;
        SingleCycleHistoHolder* singleRegularHisto, *singleIntegralHisto;
        //passed fit functions used to make the TF1 for fitting
        decayFunction* fitFunctions;
        decayFunction passedRegularFunction, passedIntegralFunction;
        TRandom3 rand;
        Double_t timeRunStart, timeRunEnd, leaveOutStartBinNumber, leaveOutEndBinNumber;
        Double_t* randArr;
        ChainFitValues* totalRegularFitParameters, *totalIntegralFitParameters;
        SingleElementFitValues* singleRegularFitParameters, *singleIntegralFitParameters;
        ParameterValue** paraVals;
        Int_t globalSeedChanger = 0;
        //helper functions
        void changeSeed();
        void createSingleFitFunctions();
        void createTotalFitFunctions();
        void DisplayParameterLimits();
        void DisplayTotalFunctionParameters();
        void genIntegralHistoSimulated();
        void setFunctionParametersTotal();
        void setFunctionParamersSingle();
        void setParaLimits();
        //TCanvas* test = new TCanvas("test", "test", 500, 500);
        clock_t Tclock;
};  


//constructor for generating a histogram
ElementFit::ElementFit(Int_t events, Int_t numRuns, Int_t numCycles, Double_t (*regularFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Double_t (**fitFunctions)(Double_t*, Double_t*), Int_t numElements,
                       Double_t timeRunEnd, Int_t numBins, string* elementNames,  ParameterValue** paraVals, Int_t singleHistoChoice, Int_t rebinChoice, Int_t rebinDifference, Double_t leaveOutStartBinNumber, Double_t leaveOutEndBinNumber)
{
    //setting variables
    this->events = events;
    this->numRuns = numRuns;
    this->numCycles = numCycles;
    this->elementNames = elementNames;
    this->fitFunctions = fitFunctions;
    this->numElements = numElements;
    this->numParameters = numElements*2;
    this->numBins = numBins;
    this->passedRegularFunction = regularFunc;
    this->passedIntegralFunction = integralFunc;
    this->timeRunEnd = timeRunEnd;
    this->paraVals = paraVals;
    this->singleHistoChoice = singleHistoChoice;
    this->rebinChoice = rebinChoice;
    this->rebinDifference = rebinDifference;
    this->leaveOutStartBinNumber = leaveOutStartBinNumber;
    this->leaveOutEndBinNumber = leaveOutEndBinNumber;

    //dynamically allocating needed variables and arrays
    //dynamic array
    singleRegularFitParameters = new SingleElementFitValues(numElements);
    singleIntegralFitParameters = new SingleElementFitValues(numElements);
    totalRegularFitParameters = new ChainFitValues(numElements);
    totalIntegralFitParameters = new ChainFitValues(numElements);
    binSizeArr = new Int_t [numCycles];
    int rebinSize = numBins;
    singleFitFunctions = new TF1* [numElements*2];
    if(rebinChoice == 1)
    {
        for(int i = 0; i < numCycles; i++)
        {
            binSizeArr[i] = rebinSize;
            rebinSize = rebinSize + rebinDifference;
        }
    }else{
        for(int i = 0; i < numCycles; i++)
        {
            binSizeArr[i] = rebinSize;
        }
    }
    randArr = new Double_t [numElements];
    timeRunStart = 0;
    Tclock = clock();
    //generates all the histo objects and 
    genAndFillHistos();
    //removes restrictions on fitting function calls and itterations
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000000);
}

ElementFit::~ElementFit()
{
    delete singleRegularFitParameters;
    delete singleIntegralFitParameters;
    delete totalRegularFitParameters;
    delete totalIntegralFitParameters;
    for(int i = 0; i < numElements*2; i++)
    {
        delete singleFitFunctions[i];
    }
    delete [] singleFitFunctions;
    delete [] randArr;
    delete regularFunction;
    delete integralFunction;
    delete regularHisto;
    delete integralHisto;
    delete singleRegularHisto;
    delete singleIntegralHisto;
    delete [] binSizeArr;
}

//used for changing seed between runs
void ElementFit::changeSeed()
{
    Double_t seederSeed;
    Double_t randNumber;

    randNumber = rand.Uniform();
    seederSeed = Tclock + globalSeedChanger + randNumber;

    rand.SetSeed(Tclock + globalSeedChanger + seederSeed);
    globalSeedChanger++;
}

//dynamically allocates the functions for the single elements
void ElementFit::createSingleFitFunctions(Double_t startFitTime, Double_t endFitTime)
{
    for(int i = 0; i < numElements; i++)
    {
        singleFitFunctions[(i*2)] = new TF1((elementNames[i] + "RegSingFunc").c_str(), fitFunctions[(i*2)], 0., timeRunEnd, (i+1)*2);
        singleFitFunctions[(i*2)+1] = new TF1((elementNames[i] + "InteSingFunc").c_str(), fitFunctions[(i*2)+1], 0., timeRunEnd, (i+1)*2);
    }
}

//creates the holding objects for the histograms
void ElementFit::createHistoHolders()
{
    string histoName;
    histoName = "Total Regular Histo";
    regularHisto = new CycleHistoHolder(numCycles, numRuns, histoName, binSizeArr, timeRunEnd);
    histoName = "Total Integral Histo";
    integralHisto = new CycleHistoHolder(numCycles, numRuns, histoName, binSizeArr, timeRunEnd);
    histoName = "Single Regular Histo";
    singleRegularHisto = new SingleCycleHistoHolder(numCycles, numElements, numRuns, histoName, binSizeArr, timeRunEnd, elementNames);
    histoName = "Single Integral Histo";
    singleIntegralHisto = new SingleCycleHistoHolder(numCycles, numElements, numRuns, histoName, binSizeArr, timeRunEnd, elementNames);
}

//dynamically allocates the total fit functions
void ElementFit::createTotalFitFunctions(Double_t startFitTime, Double_t endFitTime)
{
    regularFunction = new TF1("TotalregularFunction", passedRegularFunction, 0., timeRunEnd, numParameters);
    integralFunction = new TF1("TotalIntegralFunction", passedIntegralFunction, 0., timeRunEnd, numParameters);
}

//USED FOR TROUBLESHOOTING displays integral histogram based on canvas passed in
void ElementFit::displayIntegralHisto(TCanvas* can)
{
    can->cd();
    integralHisto->GetAHisto(0, 0)->Draw();
}

//USED FOR TROUBLESHOOTING
void ElementFit::displayParameters()
{  
    cout << endl << "Regular FIT PARAMETERS/ERRORS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ": " << "\tHalfLife: " << totalRegularFitParameters->GetAnHalfLife(i) << "s" << endl << 
        "\tError: " << totalRegularFitParameters->GetAnHalfLifeError(i) << "s" << endl;
        cout << "\tN0: " << totalRegularFitParameters->GetAnN0(i) << "" << endl << 
        "\tError: " << totalRegularFitParameters->GetAnN0Error(i) << "" << endl << endl;
    }
    cout << "INTEGRAL FIT PARAMETERS/ERRORS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ": " << "\tHalfLife: " << totalIntegralFitParameters->GetAnHalfLife(i) << "s" << endl << 
        "\tError: " << totalIntegralFitParameters->GetAnHalfLifeError(i) << "s" << endl;
        cout << "\tN0: " << totalIntegralFitParameters->GetAnN0(i) << "" << endl << 
        "\tError: " << totalIntegralFitParameters->GetAnN0Error(i) << "" << endl << endl;
    }
    cout << "REGULAR SINGLE FIT PARAMETERS/ERRORS" << endl << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << " Single Fit Values" << endl;
        for(int k = 0; k < i+1; k++)
        {
            cout << elementNames[k] << ": \tHalf Life: " << singleRegularFitParameters->GetAnHalfLife(i, k) << "s" << endl << 
            "\tHalf Life Error: " << singleRegularFitParameters->GetAnHalfLifeError(i,k) << "s" << endl;
            cout << "\tN0: " << singleRegularFitParameters->GetAnN0(i, k) << "" << endl << 
            "\tN0 Error: " << singleRegularFitParameters->GetAnN0Error(i,k) << "" << endl;
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

//USED FOR TROUBLESHOOTING displays Regular histogram
void ElementFit::displayRegularHisto(TCanvas* can)
{
    can->cd();
    regularHisto->GetAHisto(0, 0)->Draw();
}

//displays the single element histograms
void ElementFit::displaySingleHistos(TCanvas** can)
{
    for(int i = 0; i < numElements; i++)
    {
        can[(i*2)]->cd();
        singleRegularHisto->GetAHisto(0, 0, i)->Draw();
        can[(i*2)+1]->cd();
        singleIntegralHisto->GetAHisto(0, 0, i)->Draw();
    }
}

//used to display the individual histograms
void ElementFit::DrawIndividualHistos(CycleCanvasHolder* regularTotalCanvases, CycleCanvasHolder* integralTotalCanvases, SingleCycleCanvasHolder* singleRegularCanvases, SingleCycleCanvasHolder* singleIntegralCanvases,
                                      Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex)
{
    upperRunIndex = upperRunIndex + 1;
    upperCycleIndex = upperCycleIndex + 1;

    for(int cycleIndex = lowerCycleIndex; cycleIndex < upperCycleIndex; cycleIndex++)
    {
        for(int runIndex = lowerRunIndex; runIndex < upperRunIndex; runIndex++)
        {
            regularTotalCanvases->GetACanvas(cycleIndex, runIndex)->cd();
            regularHisto->GetAHisto(cycleIndex, runIndex)->Draw();

            integralTotalCanvases->GetACanvas(cycleIndex, runIndex)->cd();
            integralHisto->GetAHisto(cycleIndex, runIndex)->Draw();

            for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
            {
                singleRegularCanvases->GetACanvas(cycleIndex, runIndex, elementIndex)->cd();
                singleRegularHisto->GetAHisto(cycleIndex, runIndex, elementIndex)->Draw();

                singleIntegralCanvases->GetACanvas(cycleIndex, runIndex, elementIndex)->cd();
                singleIntegralHisto->GetAHisto(cycleIndex, runIndex, elementIndex)->Draw();
            }
        }
    }
}

//fits all histos
void ElementFit::fitHistos(Int_t cycleIndex, Int_t runIndex)
{
    createSingleFitFunctions();
    createTotalFitFunctions();
    setFunctionParametersTotal();
    setFunctionParamersSingle();
    setParaLimits();

    fitRegularHisto(cycleIndex, runIndex, leaveOutStartBinNumber, leaveOutEndBinNumber);
    fitIntegralHisto(cycleIndex, runIndex, leaveOutStartBinNumber, leaveOutEndBinNumber);
    fitSingleHistos(cycleIndex, runIndex, leaveOutStartBinNumber, leaveOutEndBinNumber);
}

//fits data of the one histogram generated
void ElementFit::fitDataGenOnce(Int_t cycleIndex, Int_t runIndex)
{
    createSingleFitFunctions();
    createTotalFitFunctions();
    setFunctionParametersTotal();
    setFunctionParamersSingle();
    setParaLimits();

    fitRegularHisto(cycleIndex, runIndex, leaveOutStartBinNumber, leaveOutEndBinNumber);
    fitIntegralHisto(cycleIndex, runIndex, leaveOutStartBinNumber, leaveOutEndBinNumber);
    fitSingleHistos(cycleIndex, runIndex, leaveOutStartBinNumber, leaveOutEndBinNumber);
}

//fits integral histogram with log likelihood method, stores fitted value in their respective arrays
void ElementFit::fitIntegralHisto(Int_t cycleIndex, Int_t runIndex, Double_t leaveOutStartBinNumber, Double_t leaveOutEndBinNumber)
{
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;
    Double_t startFit;
    Double_t endFit;
    Doublt_t binWidth = ;
    TH1D* tempHisto;

    cout << "FITTING TOTAL INTEGRAL CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;

    tempHisto = integralHisto->GetAHisto(cycleIndex, runIndex);
    tempHisto->Fit(integralFunction, "L", "", timeRunStart, timeRunEnd);
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

//fits Regular histogram with log likelihood method, stores values in their respective arrays
void ElementFit::fitRegularHisto(Int_t cycleIndex, Int_t runIndex, Double_t leaveOutStartBinNumber, Double_t leaveOutEndBinNumber)
{
    //dynamic array
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;
    TH1D* tempHisto;

    cout << "FITTING TOTAL REGULAR CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;

    tempHisto = regularHisto->GetAHisto(cycleIndex, runIndex);
    tempHisto->Fit(this->regularFunction, "L", "", timeRunStart, timeRunEnd);
    for(int i = 0; i < numElements; i++)
    {   
        valueN0 = regularFunction->GetParameter((i*2));
        errorN0 = regularFunction->GetParError((i*2));
        valueDecayConst = regularFunction->GetParameter((i*2)+1);
        errorDecayConst = regularFunction->GetParError((i*2)+1);
        valueHalfLife = log(2)/(valueDecayConst);
        errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));
        //used to get either the N0 value or half life value
        totalRegularFitParameters->SetAnN0(i, valueN0);
        totalRegularFitParameters->SetAnN0Error(i, errorN0);
        totalRegularFitParameters->SetAnHalfLife(i, valueHalfLife);
        totalRegularFitParameters->SetAnHalfLifeError(i, errorHalfLife);
    }
}

//fits all the single element histos and stores them
void ElementFit::fitSingleHistos(Int_t cycleIndex, Int_t runIndex, Double_t leaveOutStartBinNumber, Double_t leaveOutEndBinNumber)
{
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;
    TH1D* tempSingleRegularHisto, *tempSingleIntegralHisto;

    //loop for fitting and storing value for elements
    for(int i = 0; i < numElements; i++)
    {
        tempSingleRegularHisto = singleRegularHisto->GetAHisto(cycleIndex, runIndex, i);
        tempSingleIntegralHisto = singleIntegralHisto->GetAHisto(cycleIndex, runIndex, i);

        cout << "FITTING SINGLE REGULAR " << elementNames[i] << " CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;
        tempSingleRegularHisto->Fit(singleFitFunctions[(i*2)], "L", "", timeRunStart, timeRunEnd);
        cout << "FITTING SINGLE INTEGRAL " << elementNames[i] << " CYCLE: " << cycleIndex << " RUN: " << runIndex << endl;
        tempSingleIntegralHisto->Fit(singleFitFunctions[(i*2)+1], "L", "", timeRunStart, timeRunEnd);
        //loop for extracting the error and value for each element(you might need to draw out the array structure to understand what is happening)
        for(int k = 0; k < (i+1); k++)
        {
            //storing values for regular part of function
            valueN0 = singleFitFunctions[(i*2)]->GetParameter((k*2));
            errorN0 = singleFitFunctions[(i*2)]->GetParError((k*2));
            valueDecayConst = singleFitFunctions[(i*2)]->GetParameter((k*2)+1);
            errorDecayConst = singleFitFunctions[(i*2)]->GetParError((k*2)+1);
            valueHalfLife = log(2)/(valueDecayConst);
            errorHalfLife = ((log(2)/valueDecayConst)*(errorDecayConst/valueDecayConst));

            singleRegularFitParameters->SetAnN0(i, k, valueN0);
            singleRegularFitParameters->SetAnN0Error(i, k, errorN0);
            singleRegularFitParameters->SetAnHalfLife(i, k, valueHalfLife);
            singleRegularFitParameters->SetAnHalfLifeError(i, k, errorHalfLife);

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

void ElementFit::genAndFillHistos()
{
    createHistoHolders();
    genRandomAlternate();
    genIntegralHistoSimulated();
}

//generates the integral histogram from the data in the Regular histogram
void ElementFit::genIntegralHisto()
{
    //resets histogram so it can be used between runs
    TH1D* tempIntegralHisto, *tempRegularHisto;
    Int_t tempNumBins;

    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        tempNumBins = binSizeArr[cycleIndex];
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            tempIntegralHisto = integralHisto->GetAHisto(cycleIndex, runIndex);
            tempRegularHisto = regularHisto->GetAHisto(cycleIndex, runIndex);

            tempIntegralHisto->SetBinContent(1, tempRegularHisto->GetBinContent(1));
            for(int i = 2; i<= tempNumBins; i++)
            {
                tempIntegralHisto->SetBinContent(i, tempRegularHisto->GetBinContent(i) + tempIntegralHisto->GetBinContent(i-1));
            }
        }
    }
}

//generates integral histos for simulated histos
void ElementFit::genIntegralHistoSimulated()
{
    genIntegralHisto();
    genIntegralSingleHistos();
}

//generates integral for the single element hisograms
void ElementFit::genIntegralSingleHistos()
{
    TH1D* tempIntegralHisto, *tempRegularHisto;
    Int_t tempNumBins;

    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        tempNumBins = binSizeArr[cycleIndex];
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            for(int i = 0; i < numElements; i++)
            {   
                tempIntegralHisto = singleIntegralHisto->GetAHisto(cycleIndex, runIndex, i);
                tempRegularHisto = singleRegularHisto->GetAHisto(cycleIndex, runIndex, i);

                tempIntegralHisto->SetBinContent(1, tempRegularHisto->GetBinContent(1));
                for(int k = 2; k <= tempNumBins; k++)
                {
                    tempIntegralHisto->SetBinContent(k, tempIntegralHisto->GetBinContent(k-1) + tempRegularHisto->GetBinContent(k));
                }
            }
        }
    }
}

//Correct generation, Still do not know how method works yet
void ElementFit::genRandomAlternate()
{
    TH1D* tempHisto;
    TH1D** singleTempHisto;
    singleTempHisto = new TH1D* [numElements];
    Double_t hold = 0;
    Double_t stack = 0;
    //case for generating single histogram
    if(singleHistoChoice == 1 && rebinChoice == 2)
    {   //generate the single histogram
    cout << "case 1" << endl << endl;
        TH1D* originalHisto;
        TH1D** originalSingleHisto;
        originalSingleHisto = new TH1D* [numElements];

        tempHisto = regularHisto->GetAHisto(0, 0);
        for(int i = 0; i < numElements; i++)
        {
            singleTempHisto[i] = singleRegularHisto->GetAHisto(0, 0, i);
        } 
        changeSeed();
        for(int i = 0; i < events; i++)
        {
            for(int j = 0; j < numElements; j++)
            {
                hold = rand.Uniform();
                randArr[j] = (-TMath::Log(hold)) / (paraVals[j]->getValueDecayConst());
            }
            
            for(int k = 0; k < numElements; k++)
            {
                stack += randArr[k];
                singleTempHisto[k]->Fill(stack);
                tempHisto->Fill(stack);
            }
            stack = 0.0f; 
        }

        //copy generated histogram to all other histos
        originalHisto = regularHisto->GetAHisto(0, 0);
        for(int i = 0; i < numElements; i++)
        {
            originalSingleHisto[i] = singleRegularHisto->GetAHisto(0, 0, i);
        }
        for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
        {
            for(int runIndex = 0; runIndex < numRuns; runIndex++)
            {
                regularHisto->SetAHisto(cycleIndex, runIndex, originalHisto);
                for(int i = 0; i < numElements; i++)
                {
                    singleRegularHisto->SetAHisto(cycleIndex, runIndex, i, originalSingleHisto[i]);
                }
            }
        }

        delete [] originalSingleHisto;
    //case for multiple histogram generation
    }else if(singleHistoChoice == 2 && rebinChoice == 2)
    {
        cout << "case 2" << endl << endl;
        for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
        {
            for(int runIndex = 0; runIndex < numRuns; runIndex++)
            {
                tempHisto = regularHisto->GetAHisto(cycleIndex, runIndex);
                for(int i = 0; i < numElements; i++)
                {
                    singleTempHisto[i] = singleRegularHisto->GetAHisto(cycleIndex, runIndex, i);
                } 
                changeSeed();
                for(int i = 0; i < events; i++)
                {
                    for(int j = 0; j < numElements; j++)
                    {
                        hold = rand.Uniform();
                        randArr[j] = (-TMath::Log(hold)) / (paraVals[j]->getValueDecayConst());
                    }
                
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
    }else if(rebinChoice == 1)
    {
        cout << "case 3" << endl << endl;
        //numbers generated same but data fed in differently. We want all the events in each cycle to be identical so we can see the effects of rebinning, therefore data generated in run index 0 of cycle 1 must be identical to run index 0 of cycle 2
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            //change seed for generation between each run
            changeSeed();
            //data generated as normal
            for(int i = 0; i < events; i++)
            {
                for(int j = 0; j < numElements; j++)
                {
                    hold = rand.Uniform();
                    randArr[j] = (-TMath::Log(hold)) / (paraVals[j]->getValueDecayConst());
                }
                for(int k = 0; k < numElements; k++)
                {
                    stack += randArr[k];
                    //must itterate over cycles because and fill
                    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
                    {
                        tempHisto = regularHisto->GetAHisto(cycleIndex, runIndex);
                        singleTempHisto[k] = singleRegularHisto->GetAHisto(cycleIndex, runIndex, k);
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

//set parameters for base values for the total fit functions
void ElementFit::setFunctionParametersTotal()
{
    //DisplayTotalFunctionParameters();
    for(int i = 0; i < numElements; i++)
    {
        //have to have the if statments to account for the fact that the value could be a range
        if(paraVals[i]->getIsValueInitValue())
        {
            integralFunction->SetParameter((i*2), paraVals[i]->getInitValue());
            regularFunction->SetParameter((i*2), paraVals[i]->getInitValue());
        }else{
            integralFunction->SetParameter((i*2), paraVals[i]->getRangeAverageInitValue());
            regularFunction->SetParameter((i*2), paraVals[i]->getRangeAverageInitValue());
        }
        if(paraVals[i]->getIsValueDecayConst())
        {
            integralFunction->SetParameter((i*2)+1, paraVals[i]->getValueDecayConst());
            regularFunction->SetParameter((i*2)+1, paraVals[i]->getValueDecayConst());
        }else{
            integralFunction->SetParameter((i*2)+1, paraVals[i]->getRangeAverageDecayConst());
            regularFunction->SetParameter((i*2)+1, paraVals[i]->getRangeAverageDecayConst());
        }
    }
}

//displays initial parameters set for fit functions
void ElementFit::DisplayTotalFunctionParameters()
{
    cout << "INITIAL FIT FUNCTION PARAMETERS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ":\tInit value: " << paraVals[i]->getInitValue();
        cout << "\tDecay constant: " << paraVals[i]->getValueDecayConst();
        cout << endl;
    }
    cout << endl;
}

//set parameters to base values
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

//sets parameter limits so fitting knows about where to fit
void ElementFit::setParaLimits()
{
    //DisplayParameterLimits();
    //sets limits for N0 of the total function
    for(int i = 0; i < numElements; i++)
    {
        if(paraVals[i]->getIsValueInitValue())
        {
            regularFunction->SetParLimits((i*2), 0., events*2);
            integralFunction->SetParLimits((i*2), 0., events*2);
        }else{
            regularFunction->SetParLimits((i*2), paraVals[i]->getLowerRangeInitValue(), paraVals[i]->getUpperRangeInitValue());
            integralFunction->SetParLimits((i*2), paraVals[i]->getLowerRangeInitValue(), paraVals[i]->getUpperRangeInitValue());
        }
    }
    //sets limits for lambda value of total function
    for(int i = 0; i < numElements; i++)
    {
        if(paraVals[i]->getIsValueDecayConst())
        {
            regularFunction->SetParLimits((i*2)+1, (paraVals[i]->getValueDecayConst() * .001), (paraVals[i]->getValueDecayConst() * 100.));
            integralFunction->SetParLimits((i*2)+1, (paraVals[i]->getValueDecayConst() * .001), (paraVals[i]->getValueDecayConst() * 100.));
        }else{
            regularFunction->SetParLimits((i*2)+1, (paraVals[i]->getLowerRangeDecayConst()), (paraVals[i]->getUpperRangeDecayConst()));
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
                (singleFitFunctions[(i*2)])->SetParLimits((k*2), 0., events*2);
                (singleFitFunctions[(i*2)+1])->SetParLimits((k*2), 0., events*2);
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
                (singleFitFunctions[(i*2)])->SetParLimits((k*2)+1, (paraVals[k]->getValueDecayConst() * .001), (paraVals[k]->getValueDecayConst() * 100.));
                (singleFitFunctions[(i*2)+1])->SetParLimits((k*2)+1, (paraVals[k]->getValueDecayConst() * .001), (paraVals[k]->getValueDecayConst() * 100.));
            }else{
                (singleFitFunctions[(i*2)])->SetParLimits((k*2)+1, (paraVals[i]->getLowerRangeDecayConst()), (paraVals[i]->getUpperRangeDecayConst()));
                (singleFitFunctions[(i*2)+1])->SetParLimits((k*2)+1, (paraVals[i]->getLowerRangeDecayConst()), (paraVals[i]->getUpperRangeDecayConst()));
            }
        }
    }
}

void ElementFit::DisplayParameterLimits()
{
    cout << "PARAMETER FIT RANGE" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ":" << endl;
        cout << "\tInitial Lower Range: 0" << endl;
        cout << "\tInitial Lower Range: " << events * 2 << endl;
        cout << "\tDecay Constant Lower Range: " << paraVals[i]->getValueDecayConst() * .001 << endl;
        cout << "\tDecay Constant Upper Range: " << paraVals[i]->getValueDecayConst() * 100. << endl;
        cout << endl;
    }
}

#endif