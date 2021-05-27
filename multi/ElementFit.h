#ifndef ELEMENTFIT_H
#define ELEMENTFIT_H

#include <ctime>
#include <iostream>
#include <string>
#include <math.h>

#include "TF1.h"
#include "TH1.h"
#include "time.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TList.h"
#include "ParameterValue.h"
#include "FitParameterStore.h"
#include "ChainFitValues.h"
#include "FitValues.h"
#include "SingleElementFitValues.h"

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
        TH1D* getIntegralHisto(){return integralHisto;}
        TH1D* getRegularHisto(){return regularHisto;}
        Int_t getNumEvents(){return events;}
        Int_t getNumElements(){return numElements;}
        Int_t getNumBins(){return numBins;}
        Int_t getTimeRunEnd(){return timeRunEnd;}
        ChainFitValues* getRegularFitParameters(){return totalRegularFitParameters;}
        ChainFitValues* getIntegralFitParameters(){return totalIntegralFitParameters;}
        SingleElementFitValues* getSingleRegularFitParameters(){return singleRegularFitParameters;}
        SingleElementFitValues* getSingleIntegralFitParameters(){return singleIntegralFitParameters;}
        //setter function
        void setNumBins(Int_t numBins);
        void setNumEvents(Int_t events){this->events = events;} 
        void setTimeRunEnd(Double_t timeRunEnd){this->timeRunEnd = timeRunEnd;}
        void setTimeRunStart(Double_t timeRunStart){this->timeRunStart = timeRunStart;}
        //public functions
        void createHistoFile(string name);
        void displayIntegralHisto(TCanvas* can);
        void displayRegularHisto(TCanvas* can);
        void displaySingleHistos(TCanvas** can);
        void displayParameters();
        void fitData();
        void fitDataGenOnce();
        void fitDataLoadedHisto();
        void genIntegralHisto();
        void genIntegralSingleHistos();
        void genRandomAlternate();
        void writeToFile();
        void fitRegularHisto();
        void fitIntegralHisto();
    private:
        //private variables
        string* elementNames;
        Int_t events, numElements, numBins, numParameters, doubleNumElements;
        TF1* integralFunction, *regularFunction;
        TF1** singleFitFunctions;
        TH1D** singleElementRegularHistos, **singleElementIntegralHistos;
        TH1D* regularHisto, *integralHisto;
        TList* histoList;
        TFile* histoFile;
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
        void createSingleElementRegularHistos();
        void createSingleElementIntegralHistos();
        void createTotalFunctions();
        void createTotalHistos();
        void fitSingleHistos();
        void genIntegralHistoSimulated();
        void setFunctionParametersTotal();
        void setFunctionParamersSingle();
        void setParaLimits();
        //TCanvas* test = new TCanvas("test", "test", 500, 500);
        clock_t Tclock;
        TCanvas *CSCycleResultCanvas;
        TCanvas *BACycleResultCanvas;
        TCanvas *LACycleResultCanvas;
        TCanvas *CSCycleResultCanvasInt;
        TCanvas *BACycleResultCanvasInt;
        TCanvas *LACycleResultCanvasInt;
        TCanvas** testCanArr;
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
    /*
    CSCycleResultCanvas = new TCanvas("CSCycleResultCanvas", "CSCycleResultCanvas", 500, 500);
    BACycleResultCanvas = new TCanvas("BACycleResultCanvas", "BACycleResultCanvas", 500, 500);
    LACycleResultCanvas = new TCanvas("LACycleResultCanvas", "LACycleResultCanvas", 500, 500);
    CSCycleResultCanvasInt = new TCanvas("CSCycleResultCanvasInt", "CSCycleResultCanvasInt", 500, 500);
    BACycleResultCanvasInt = new TCanvas("BACycleResultCanvasInt", "BACycleResultCanvasInt", 500, 500);
    LACycleResultCanvasInt = new TCanvas("LACycleResultCanvasInt", "LACycleResultCanvasInt", 500, 500);
    testCanArr = new TCanvas* [6];
    testCanArr[0] = CSCycleResultCanvas;
    testCanArr[1] = BACycleResultCanvas;
    testCanArr[2] = LACycleResultCanvas;
    testCanArr[3] = CSCycleResultCanvasInt;
    testCanArr[4] = BACycleResultCanvasInt;
    testCanArr[5] = LACycleResultCanvasInt;
    */

    //dynamically allocating needed variables and arrays
    histoList = new TList();
    //dynamic array
    singleRegularFitParameters = new SingleElementFitValues(numElements);
    singleIntegralFitParameters = new SingleElementFitValues(numElements);
    totalRegularFitParameters = new ChainFitValues(numElements);
    totalIntegralFitParameters = new ChainFitValues(numElements);
    singleElementRegularHistos = new TH1D* [numElements];
    singleElementIntegralHistos = new TH1D* [numElements];
    randArr = new Double_t [numElements];
    timeRunStart = 0;
    Tclock = clock();
    //dynamically allocates single element histos
    createSingleElementRegularHistos();
    //dynamically allocates histograms for the total histogram
    createTotalHistos();
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
                       ParameterValue** paraVals_p, TH1D* loadedHisto_p, string* elementNames_p)
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
    delete [] singleElementRegularHistos;
    delete [] singleElementIntegralHistos;
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
    delete histoList;
    delete regularFunction;
    delete integralFunction;
}

//used for changing seed between runs
void ElementFit::changeSeed()
{
    rand.SetSeed(Tclock + globalSeedChanger);
    globalSeedChanger++;
}

//creates file to export histograms to
void ElementFit::createHistoFile(string name)
{
    histoFile = new TFile((name).c_str(), "recreate");
}

//creates the total histograms
void ElementFit::createTotalHistos()
{
    regularHisto = new TH1D("TotalRegularHisto", "TotalRegularHisto", numBins, 0., timeRunEnd);
    regularHisto->GetXaxis()->SetTitle("Time (S)");
    regularHisto->GetYaxis()->SetTitle("Counts");
    integralHisto = new TH1D("TotalIntegralHisto", "TotalIntegralHisto", numBins, 0., timeRunEnd);
    integralHisto->GetXaxis()->SetTitle("Time (S)");
    integralHisto->GetYaxis()->SetTitle("Counts");
    histoList->Add(regularHisto);
    histoList->Add(integralHisto);
}

//creates the histograms for the single elements 
void ElementFit::createSingleElementRegularHistos()
{
    for(int i = 0; i < numElements; i++)
    {   
        singleElementRegularHistos[i] = new TH1D(((elementNames[i] + "SingleRegularHisto").c_str()), ((elementNames[i] + "SingleRegularHisto").c_str()), numBins, 0., timeRunEnd);
        singleElementRegularHistos[i]->GetXaxis()->SetTitle("Time(s)");
        singleElementRegularHistos[i]->GetYaxis()->SetTitle("Number of Events");
        histoList->Add(singleElementRegularHistos[i]);
        singleElementIntegralHistos[i] = new TH1D(((elementNames[i] + "SingleIntegralHisto").c_str()), ((elementNames[i] + "SingleIntegralHisto").c_str()), numBins, 0., timeRunEnd);
        singleElementIntegralHistos[i]->GetXaxis()->SetTitle("Time(s)");
        singleElementIntegralHistos[i]->GetYaxis()->SetTitle("Number of Events");
        histoList->Add(singleElementIntegralHistos[i]);
    }
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
    integralHisto->Draw();
}

//USED FOR TROUBLESHOOTING
void ElementFit::displayParameters()
{  
    cout << endl << "Regular FIT PARAMETERS/ERRORS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ": " << "\tHalfLife: " << totalRegularFitParameters->GetAnHalfLife(i) << "s" << endl << 
        "\tError: " << totalRegularFitParameters->GetAnHalfLifeError(i) << "s" << endl;
    }
    cout << endl << "INTEGRAL FIT PARAMETERS/ERRORS" << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << ": " << "\tHalfLife: " << totalIntegralFitParameters->GetAnHalfLife(i) << "s" << endl << 
        "\tError: " << totalIntegralFitParameters->GetAnHalfLifeError(i) << "s" << endl;
    }
    cout << endl << "REGULAR SINGLE FIT PARAMETERS/ERRORS" << endl << endl;
    for(int i = 0; i < numElements; i++)
    {
        cout << elementNames[i] << " Single Fit Values" << endl;
        for(int k = 0; k < i+1; k++)
        {
            cout << elementNames[k] << ": \tHalf Life: " << singleRegularFitParameters->GetAnHalfLife(i, k) << "s" << endl << 
            "\tHalf Life Error: " << singleRegularFitParameters->GetAnHalfLifeError(i,k) << "s" << endl;
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
            cout << endl;
        }
    }
}

//USED FOR TROUBLESHOOTING displays Regular histogram
void ElementFit::displayRegularHisto(TCanvas* can)
{
    can->cd();
    regularHisto->Draw();
}

//displays the single element histograms
void ElementFit::displaySingleHistos(TCanvas** can)
{
    for(int i = 0; i < numElements; i++)
    {
        can[(i*2)]->cd();
        singleElementRegularHistos[i]->Draw();
        can[(i*2)+1]->cd();
        singleElementIntegralHistos[i]->Draw();
    }
}

//generates data for the integral and Regular histograms and then fits them using helper functions
void ElementFit::fitData()
{
    genRandomAlternate();
    genIntegralHistoSimulated();
    fitRegularHisto();
    fitIntegralHisto();
    fitSingleHistos();
}

//fits data of the one histogram generated
void ElementFit::fitDataGenOnce()
{
    fitRegularHisto();
    fitIntegralHisto();
    fitSingleHistos();
}

//fit data for loaded in histogram
void ElementFit::fitDataLoadedHisto()
{
    genIntegralHisto();
    fitRegularHisto();
    fitIntegralHisto();
}

//fits integral histogram with log likelihood method, stores fitted value in their respective arrays
void ElementFit::fitIntegralHisto()
{
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;

    integralHisto->Fit(integralFunction, "L", "", timeRunStart, timeRunEnd);
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
void ElementFit::fitRegularHisto()
{
    //dynamic array
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;

    regularHisto->Fit(this->regularFunction, "L", "", timeRunStart, timeRunEnd);
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
void ElementFit::fitSingleHistos()
{
    Double_t valueN0;
    Double_t errorN0;
    Double_t valueDecayConst;
    Double_t errorDecayConst;
    Double_t valueHalfLife;
    Double_t errorHalfLife;

    //loop for fitting and storing value for elements
    for(int i = 0; i < numElements; i++)
    {
        singleElementRegularHistos[i]->Fit(singleFitFunctions[(i*2)], "L", "", timeRunStart, timeRunEnd);
        singleElementIntegralHistos[i]->Fit(singleFitFunctions[(i*2)+1], "L", "", timeRunStart, timeRunEnd);
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
    integralHisto->Reset("ICES");
    integralHisto->SetBinContent(1, regularHisto->GetBinContent(1));
    for(int i = 2; i<= numBins; i++)
    {
        integralHisto->SetBinContent(i, regularHisto->GetBinContent(i) + integralHisto->GetBinContent(i-1));
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
    for(int i = 0; i < numElements; i++)
    {
        singleElementIntegralHistos[i]->Reset("ICES");
    }

    for(int i = 0; i < numElements; i++)
    {   
        singleElementIntegralHistos[i]->SetBinContent(1, singleElementRegularHistos[i]->GetBinContent(1));
        for(int k = 2; k <= events; k++)
        {
            singleElementIntegralHistos[i]->SetBinContent(k, singleElementIntegralHistos[i]->GetBinContent(k-1) + singleElementRegularHistos[i]->GetBinContent(k));
        }
    }
}

//Correct generation, Still do not know how method works yet
void ElementFit::genRandomAlternate()
{
    Double_t hold = 0;
    regularHisto->Reset("ICES");
    for(int i = 0; i < numElements; i++)
    {
        singleElementRegularHistos[i]->Reset("ICES");
    }
    changeSeed();
    for(int i = 0; i < events; i++)
    {
        for(int j = 0; j < numElements; j++)
        {
            hold = rand.Uniform();
            randArr[j] = (-TMath::Log(hold)) / (paraVals[j]->getValueDecayConst());
        }
    
        hold = 0;
        for(int k = 0; k < numElements; k++)
        {
            hold += randArr[k];
            singleElementRegularHistos[k]->Fill(hold);
            regularHisto->Fill( hold );
        }
    }
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

//changes number of bins in the histograms at runtime
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

void ElementFit::writeToFile()
{
    histoList->Write();
    histoFile->Close();
}
#endif