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
#include "RunHistoHolder.h"
#include "SingleRunHistoHolder.h"
#include "CycleHistoHolder.h"

using namespace std;

typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);

class ElementFit{
    public:
        ElementFit(Int_t events_p, Double_t (*regularFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Double_t (**fitFunctions)(Double_t*, Double_t*),
                   Int_t numElements_p, Double_t timeRunEnd_p, Int_t numBins_p, string* elementNames_p,  ParameterValue** paraVals_p);
        ElementFit(Double_t (*regularFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Int_t numElements_p, Double_t timeRunEnd_p,
                   ParameterValue** paraVals_p, TH1D* loadedHisto_p, string* elementNames_p);
        ~ElementFit();
        //getter function
        Double_t getElementParameters(int i){return paraVals[i]->getValueDecayConst();}
        Int_t getNumEvents(){return events;}
        Int_t getNumElements(){return numElements;}
        Int_t getNumBins(){return numBins;}
        Int_t getTimeRunEnd(){return timeRunEnd;}
        ChainFitValues* getRegularFitParameters(){return totalRegularFitParameters;}
        ChainFitValues* getIntegralFitParameters(){return totalIntegralFitParameters;}
        SingleElementFitValues* getSingleRegularFitParameters(){return singleRegularFitParameters;}
        SingleElementFitValues* getSingleIntegralFitParameters(){return singleIntegralFitParameters;}
        //setter function
        //void setNumBins(Int_t numBins);
        void setNumEvents(Int_t events){this->events = events;} 
        void setTimeRunEnd(Double_t timeRunEnd){this->timeRunEnd = timeRunEnd;}
        void setTimeRunStart(Double_t timeRunStart){this->timeRunStart = timeRunStart;}
        void setNumRuns(Int_t numRuns){this->numRuns = numRuns;}
        void setNumCycles(Int_t numCycles){this->numCycles = numCycles;}
        //public functions
        void createHistoHolders();
        void displayIntegralHisto(TCanvas* can);
        void displayRegularHisto(TCanvas* can);
        void displaySingleHistos(TCanvas** can);
        void displayParameters();
        
        void fitHistos(Int_t cycleIndex, Int_t runIndex);
        void fitSingleHisto(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex);
        void fitDataGenOnce();
        void fitDataLoadedHisto();
        
        void genIntegralHisto();
        void genIntegralSingleHistos();
        void genRandomAlternate();
        void genAndFillHistos();
        void fitRegularHisto(Int_t cycleIndex, Int_t runIndex);
        void fitIntegralHisto(Int_t cycleIndex, Int_t runIndex);
    private:
        //private variables
        string* elementNames;
        Int_t events, numElements, numBins, numParameters, numRuns, numCycles, singleHistoChoice;
        TF1* integralFunction, *regularFunction;
        TF1** singleFitFunctions;
        CycleHistoHolder* regularHisto, *integralHisto;
        SingleCycleHistoHolder* singleRegularHisto, *singleIntegralHisto;
        //passed fit functions used to make the TF1 for fitting
        decayFunction* fitFunctions;
        decayFunction passedRegularFunction, passedIntegralFunction;
        TRandom3 rand;
        Double_t timeRunStart, timeRunEnd;
        Double_t* randArr;
        ChainFitValues* totalRegularFitParameters;
        ChainFitValues* totalIntegralFitParameters;
        SingleElementFitValues* singleRegularFitParameters;
        SingleElementFitValues* singleIntegralFitParameters;
        ParameterValue** paraVals;
        Int_t globalSeedChanger = 0;
        //helper functions
        void changeSeed();
        TF1** createFitFunctions();
        void createTotalFunctions();
        void fitSingleHistos(Int_t cycleIndex, Int_t run);
        void genIntegralHistoSimulated();
        void setFunctionParametersTotal();
        void setFunctionParamersSingle();
        void setParaLimits();
        //TCanvas* test = new TCanvas("test", "test", 500, 500);
        clock_t Tclock;
};  


//constructor for generating a histogram
ElementFit::ElementFit(Int_t events_p, Double_t (*regularFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Double_t (**fitFunctions)(Double_t*, Double_t*),
                       Int_t numElements_p, Double_t timeRunEnd_p, Int_t numBins_p, string* elementNames_p,  ParameterValue** paraVals_p)
{
    //setting variables
    this->events = events_p;
    this->elementNames = elementNames_p;
    this->fitFunctions = fitFunctions;
    this->numElements = numElements_p;
    this->numParameters = numElements*2;
    this->numBins = numBins_p;
    this->passedRegularFunction = regularFunc;
    this->passedIntegralFunction = integralFunc;
    this->timeRunEnd = timeRunEnd_p;
    this->paraVals = paraVals_p;
    numRuns = 1;
    numCycles = 1;
    singleHistoChoice = 0;

    //dynamically allocating needed variables and arrays
    //dynamic array
    singleRegularFitParameters = new SingleElementFitValues(numElements);
    singleIntegralFitParameters = new SingleElementFitValues(numElements);
    totalRegularFitParameters = new ChainFitValues(numElements);
    totalIntegralFitParameters = new ChainFitValues(numElements);
    randArr = new Double_t [numElements];
    timeRunStart = 0;
    Tclock = clock();
    //dynamically allocates the fit functions for the single histograms
    singleFitFunctions = createFitFunctions();
    //dynamically allocates histograms for the single histograms
    createTotalFunctions();
    //setting base parameters for ALL fit functions
    setFunctionParametersTotal();
    setFunctionParamersSingle();
    //parameter limits so we get reasonable values
    setParaLimits();
}


//constructor for loading in a histogram(NOT MODIFIED WITH NEW STRUCTURE IN MIND)
/*
ElementFit::ElementFit(Double_t (*regularFunc)(Double_t*, Double_t*), Double_t (*integralFunc)(Double_t*, Double_t*), Int_t numElements_p, Double_t timeRunEnd_p,
                       ParameterValue** paraVals_p, TH1D* loadedHisto_p, string* elementNames_p, Int_t numFits)
{
    this->passedRegularFunction = regularFunc;
    this->passedIntegralFunction = integralFunc;
    this->timeRunEnd = timeRunEnd_p;
    this->paraVals = paraVals_p;
    this->numElements = numElements_p;
    this->regularHisto = loadedHisto_p;
    this->elementNames = elementNames_p;
    numBins = regularHisto->GetNbinsX();

    totalFitParameters = new FitParameterStore(numElements);
    randArr = new Double_t [numElements];
    Tclock = clock();
    timeRunStart = 0;
    integralHisto = new TH1D((elementNames[numElements-1] + "IntegralHisto").c_str(), (elementNames[numElements-1] + "IntegralHisto").c_str(), numBins, 0., timeRunEnd);

    //dynamically allocates histograms for the single histograms
    createTotalFunctions();
    //setting base parameters for ALL fit functions
    setFunctionParametersTotal();
    //parameter limits so we get reasonable values
    setParaLimits();
}
*/

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
}

//used for changing seed between runs
void ElementFit::changeSeed()
{
    rand.SetSeed(Tclock + globalSeedChanger);
    globalSeedChanger++;
}

//dynamically allocates the functions for the single elements
TF1** ElementFit::createFitFunctions()
{
    TF1** tempFunctionHolder = new TF1* [numElements*2];
    for(int i = 0; i < numElements; i++)
    {
        tempFunctionHolder[(i*2)] = new TF1((elementNames[i] + "RegSingFunc").c_str(), fitFunctions[(i*2)], 0., timeRunEnd, (i+1)*2);
        tempFunctionHolder[(i*2)+1] = new TF1((elementNames[i] + "InteSingFunc").c_str(), fitFunctions[(i*2)+1], 0., timeRunEnd, (i+1)*2);
    }
    return tempFunctionHolder;
}

//creates the holding objects for the histograms
void ElementFit::createHistoHolders()
{
    string histoName;
    histoName = "Total Regular Histo";
    regularHisto = new CycleHistoHolder(numCycles, numRuns, histoName, numBins, timeRunEnd);
    histoName = "Total Integral Histo";
    integralHisto = new CycleHistoHolder(numCycles, numRuns, histoName, numBins, timeRunEnd);
    histoName = "Single Regular Histo";
    singleRegularHisto = new SingleCycleHistoHolder(numCycles, numElements, numRuns, histoName, numBins, timeRunEnd, elementNames);
    histoName = "Single Integral Histo";
    singleIntegralHisto = new SingleCycleHistoHolder(numCycles, numElements, numRuns, histoName, numBins, timeRunEnd, elementNames);
}

//dynamically allocates the total fit functions
void ElementFit::createTotalFunctions()
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
        "\tError: " << totalRegularFitParameters->GetAnN0Error(i) << "" << endl;
    }
    cout << endl << "INTEGRAL FIT PARAMETERS/ERRORS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ": " << "\tHalfLife: " << totalIntegralFitParameters->GetAnHalfLife(i) << "s" << endl << 
        "\tError: " << totalIntegralFitParameters->GetAnHalfLifeError(i) << "s" << endl;
        cout << "\tN0: " << totalIntegralFitParameters->GetAnN0(i) << "" << endl << 
        "\tError: " << totalIntegralFitParameters->GetAnN0Error(i) << "" << endl;
    }
    cout << endl << "REGULAR SINGLE FIT PARAMETERS/ERRORS" << endl << endl;
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

void ElementFit::genAndFillHistos()
{
    if(singleHistoChoice == 1)
    {
        createHistoHolders();
    }else{
        createHistoHolders();
        genRandomAlternate();
        genIntegralHistoSimulated();
    }
}

//fits all histos
void ElementFit::fitHistos(Int_t cycleIndex, Int_t runIndex)
{
    fitRegularHisto(cycleIndex, runIndex);
    fitIntegralHisto(cycleIndex, runIndex);
    fitSingleHistos(cycleIndex, runIndex);
}

//fits data of the one histogram generated
void ElementFit::fitDataGenOnce()
{
    //fitRegularHisto();
    //fitIntegralHisto();
    //fitSingleHistos();
}

//fit data for loaded in histogram(not updated)
void ElementFit::fitDataLoadedHisto()
{
    //genIntegralHisto();
    //fitRegularHisto();
    //fitIntegralHisto();
}

//fits integral histogram with log likelihood method, stores fitted value in their respective arrays
void ElementFit::fitIntegralHisto(Int_t cycleIndex, Int_t runIndex)
{
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;
    TH1D* tempHisto;

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
void ElementFit::fitRegularHisto(Int_t cycleIndex, Int_t runIndex)
{
    //dynamic array
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;
    TH1D* tempHisto;

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
void ElementFit::fitSingleHistos(Int_t cycleIndex, Int_t runIndex)
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

        tempSingleRegularHisto->Fit(singleFitFunctions[(i*2)], "L", "", timeRunStart, timeRunEnd);
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

//generates the integral histogram from the data in the Regular histogram
void ElementFit::genIntegralHisto()
{
    //resets histogram so it can be used between runs
    TH1D* tempIntegralHisto, *tempRegularHisto;
    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            tempIntegralHisto = integralHisto->GetAHisto(cycleIndex, runIndex);
            tempRegularHisto = regularHisto->GetAHisto(cycleIndex, runIndex);

            tempIntegralHisto->SetBinContent(1, tempRegularHisto->GetBinContent(1));
            for(int i = 2; i<= numBins; i++)
            {
                tempIntegralHisto->SetBinContent(i, tempRegularHisto->GetBinContent(i) + tempIntegralHisto->GetBinContent(i-1));
            }

            ofstream myFile;
            Double_t binData;
            myFile.open("integralValues.txt");
            for(int i = 0; i < numBins; i++)
            {
                binData = tempIntegralHisto->GetBinContent(i+1);
                myFile << binData;
                myFile << "\n";
            }
            myFile.close();
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
    for(int cycleIndex = 0; cycleIndex < numCycles; cycleIndex++)
    {
        for(int runIndex = 0; runIndex < numRuns; runIndex++)
        {
            for(int i = 0; i < numElements; i++)
            {   
                tempIntegralHisto = singleIntegralHisto->GetAHisto(cycleIndex, runIndex, i);
                tempRegularHisto = singleRegularHisto->GetAHisto(cycleIndex, runIndex, i);

                tempIntegralHisto->SetBinContent(1, tempRegularHisto->GetBinContent(1));
                for(int k = 2; k <= events; k++)
                {
                    tempIntegralHisto->SetBinContent(k, tempIntegralHisto->GetBinContent(k-1) + tempRegularHisto->GetBinContent(k));
                }
            }

            ofstream myFile;
            Double_t binData;
            myFile.open("CsIntegralValues.txt");
                for(int i = 0; i < numBins; i++)
                {
                    binData = singleIntegralHisto->GetAHisto(0, 0, 0)->GetBinContent(i+1);
                    myFile << binData;
                    myFile << "\n";
                }
                myFile.close();

            myFile.open("BaIntegralValues.txt");
                for(int i = 0; i < numBins; i++)
                {
                    binData = singleIntegralHisto->GetAHisto(0, 0, 1)->GetBinContent(i+1);
                    myFile << binData;
                    myFile << "\n";
                }
                myFile.close();
            
            myFile.open("LaIntegralValues.txt");
                for(int i = 0; i < numBins; i++)
                {
                    binData = singleIntegralHisto->GetAHisto(0, 0, 2)->GetBinContent(i+1);
                    myFile << binData;
                    myFile << "\n";
                }
                myFile.close();
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
    if(singleHistoChoice == 1)
    {   //generate the single histogram
        TH1D* originalHisto;
        TH1D**originalSingleHisto;
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
        /*
        originalHisto = regularHisto->GetAHisto(0, 0);
        for(int i = 0; i < numElements; i++)
        {
            originalSingleHisto = singleRegularHisto(0, 0, i);
        }
        for(int cycleIndex = 1; cycleIndex < numCycles; cycleIndex++)
        {
            for(int runIndex = 1; runIndex < numRuns; runIndex++)
            {
                
            }
        }
        */
    }else{
    //case for multiple histogram generation
    string fileNames [3] = {"CsRegularValues.txt", "BaRegularValues.txt", "LaRegularValues.txt"};
    ofstream myFile;
    Double_t binData;
    for(int i = 0; i < numElements; i++)
    {
        myFile.open(fileNames[i]);
        myFile << "";
        myFile.close();
    }
    myFile.open("regularValues.txt");
    myFile << "";
    myFile.close();
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

                        myFile.open(fileNames[k], ios::out | ios::app);
                        myFile << stack;
                        myFile << endl;
                        myFile.close();
                        myFile.open("regularValues.txt", ios::out | ios::app);
                        myFile << stack;
                        myFile << endl;
                        myFile.close();
                    }
                    stack = 0.0f;
                }
            }
        }
    }

    delete [] singleTempHisto;
}

//set parameters for base values for the total fit functions
void ElementFit::setFunctionParametersTotal()
{
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

//changes number of bins in the histograms at runtime(not updated)
/*
void ElementFit::setNumBins(Int_t numBins)
{

    //we must specify every single bin edge so that is done here
    Double_t* binEdges = new Double_t[numBins+1];
    Double_t binWidth = (timeRunEnd/(Double_t)numBins);
    binEdges[0] = 0.;
    for(int i = 1; i < numBins; i++)
    {
        binEdges[i] = binEdges[i-1] + binWidth;
    }
    binEdges[numBins] = timeRunEnd;

    //uses rebin function in root to rebin both histograms
    regularHisto = (TH1D*)regularHisto->Rebin(numBins, (elementNames[numElements-1] + "regularHisto").c_str(), binEdges);
    integralHisto = (TH1D*)integralHisto->Rebin(numBins, (elementNames[numElements-1] + "IntegralHisto").c_str(), binEdges);
    this->numBins = numBins;
    delete [] binEdges;
}
*/

//sets parameter limits so fitting knows about where to fit
void ElementFit::setParaLimits()
{
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

#endif