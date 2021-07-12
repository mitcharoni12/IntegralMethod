#ifndef RUN_H
#define RUN_H

#include <iostream>
#include <string>
#include <iomanip>
#include <utility>
#include <vector>
#include <cstdio>
#include <sys/stat.h>

#include "TH1.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"

#include "ElementFit.h"
#include "SingleElementFitValues.h"
#include "ChainFitValues.h"
#include "SingleChainRunFitValues.h"
#include "ChainRunFitValues.h"
#include "FitOption.h"

using namespace std;

class Run{
    private:
        Int_t runs, numElements, eventDecrement;
        FitOption* fitOptions;
        ElementFit* element;
        ChainRunFitValues* batemanFitValues, *integralFitValues;
        SingleChainRunFitValues* singleBatemanFitValues, *singleIntegralFitValues;
        Double_t* eventsXAxis, *runsXAxis, *zero;
        string* elementNameStrs;
        TH3D* correlationHisto;
        TH1D** multiRunResultHistograms;
        TGraphErrors** totalBatemanGraphs, **totalIntegralGraphs, **singleBatemanGraphs, **singleIntegralGraphs;
        TCanvas* testCan;
        TCanvas** multiRunResultCanvas;
        //helper functions
        void createCorrelationHistoEventChange();
        void createCorrelationHistoNoChange();
        Double_t getMaxElement(Double_t* arr);
        Double_t getMinElement(Double_t* arr);
        void runEventChange();
        void genIntegralMeanGraphs();
        void genBatemanMeanGraphs();
    public:
        Run(ElementFit* element);
        ~Run();
        void runNoChangeGenOnce(Int_t cycleIndex, Int_t runIndex);
        void genGraphsEventChange();
        void genGraphsNoChange();
        void genGraphsNoChangeSingleElement();
        TH1D** createRunResultHistos();
        TH1D** createRunResultHistosSingleElements();
        TH1D** fillRunResultHistos(TH1D** multiRunResultHistograms);
        TH1D** fillRunResultHistosSingleElement(TH1D** multiRunResultHistogramsSingleElement);
        void displayMultiRunResultGraphs(TCanvas** canvasArray);
        void displayCorrelationHisto(TCanvas* canvas);
        void displayMultiRunResultHistos(TCanvas** canvasArray, TH1D** multiRunResultHistograms);
        void runNoChange(Int_t cycleIndex);
        //getter function
        ChainRunFitValues* getBatemanFitValues(){return batemanFitValues;}
        ChainRunFitValues* getIntegralFitValues(){return integralFitValues;}
        SingleChainRunFitValues* getSingleBatemanFitValues(){return singleBatemanFitValues;}
        SingleChainRunFitValues* getSingleIntegralFitValues(){return singleIntegralFitValues;}
        TH1D** getMultiRunResultHistos(){return multiRunResultHistograms;}
        string* getElementStringNames(){return elementNameStrs;}
        Int_t getNumRuns(){return runs;}
        //setter function
        void setNumRuns(Int_t numRuns){this->runs = numRuns;}
};

Run::Run(ElementFit* element)
{
    //variable parameter setting
    this->fitOptions = element->getFitOptions();
    this->runs = fitOptions->GetNumRuns();
    this->eventDecrement = fitOptions->GetEventDecrement();
    this->element = element;
    this->elementNameStrs = fitOptions->GetElementNames();
    numElements = fitOptions->GetNumElements();

    //dynamically allocating required arrays and root variables
    //dynamic array
    eventsXAxis = new Double_t [runs];
    runsXAxis = new Double_t [runs];
    zero = new Double_t [runs];
    for(int i = 0; i < runs; i++)
    {
        zero[i] = 0.0f;
    }
    batemanFitValues = new ChainRunFitValues(numElements, runs);
    integralFitValues = new ChainRunFitValues(numElements, runs);
    singleBatemanFitValues = new SingleChainRunFitValues(numElements, runs);
    singleIntegralFitValues = new SingleChainRunFitValues(numElements, runs);
    totalBatemanGraphs = new TGraphErrors* [numElements];
    totalIntegralGraphs = new TGraphErrors* [numElements];
    singleBatemanGraphs = new TGraphErrors* [numElements];
    singleIntegralGraphs = new TGraphErrors* [numElements];
    //testCan = new TCanvas("testCan", "testCan", 500, 500);
    for(int i = 0; i < runs; i++)
    {
        runsXAxis[i] = i+1;
    }
}

Run::~Run()
{
    delete batemanFitValues;
    delete integralFitValues;
    delete singleBatemanFitValues;
    delete singleIntegralFitValues;
    delete [] eventsXAxis;
    delete [] runsXAxis;
    delete [] totalBatemanGraphs;
    delete [] totalIntegralGraphs;
    delete [] singleBatemanGraphs;
    delete [] singleIntegralGraphs;
}

//creates histograms for the run result of the total functions
TH1D** Run::createRunResultHistos()
{
    TH1D** multiRunResultHistograms = new TH1D* [numElements*2];
    Double_t parameterValue;

    //creates the histograms
    for(int i = 0; i < numElements; i++)
    {
        parameterValue = TMath::LogE()/(element->getElementParameters(i));
        multiRunResultHistograms[(i*2)] = new TH1D((elementNameStrs[i] + " Fit Result Bateman Histo").c_str(), (elementNameStrs[i] + " Fit Result Bateman Histo").c_str(), 500, parameterValue*0, parameterValue*2.5);
        multiRunResultHistograms[(i*2)+1] = new TH1D((elementNameStrs[i] + " Fit Result Integral Histo").c_str(), (elementNameStrs[i] + " Fit Result Integral Histo").c_str(), 500, parameterValue*0, parameterValue*2.5);
    }

    return multiRunResultHistograms;
}

//creates histograms for the run result of the single element functions
TH1D** Run::createRunResultHistosSingleElements()
{
    TH1D** multiRunResultHistosSingleElement = new TH1D* [numElements*2];
    Double_t parameterValue;

    //creates the histograms
    for(int i = 0; i < numElements; i++)
    {
        parameterValue = TMath::LogE()/(element->getElementParameters(i));
        multiRunResultHistosSingleElement[(i*2)] = new TH1D((elementNameStrs[i] + " Fit Result Bateman Histo Single Element").c_str(), (elementNameStrs[i] + " Fit Result Bateman Histo Single Element").c_str(), 500, parameterValue*0, parameterValue*2.5);
        multiRunResultHistosSingleElement[(i*2)+1] = new TH1D((elementNameStrs[i] + " Fit Result Integral Histo Single Element").c_str(), (elementNameStrs[i] + " Fit Result Integral Histo Single Element").c_str(), 500, parameterValue*0, parameterValue*2.5);
    }

    return multiRunResultHistosSingleElement;
}

//display the correlation histogram
void Run::displayCorrelationHisto(TCanvas* canvas)
{
    canvas->cd();
    correlationHisto->Draw("lego");
}

//displays the graphs for the multiple runs
void Run::displayMultiRunResultGraphs(TCanvas** canvasArray)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(1);
        totalBatemanGraphs[i]->Draw();
        canvasArray[i]->cd(2);
        totalIntegralGraphs[i]->Draw();
        canvasArray[i]->cd(3);
        singleBatemanGraphs[i]->Draw();
        canvasArray[i]->cd(4);
        singleIntegralGraphs[i]->Draw();
    }
}

//displays the histograms for the multiple runs
void Run::displayMultiRunResultHistos(TCanvas** canvasArray, TH1D** multiRunResultHistograms)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(1);
        multiRunResultHistograms[(i*4)]->Draw();
        canvasArray[i]->cd(2);
        multiRunResultHistograms[(i*4)+1]->Draw();
        canvasArray[i]->cd(3);
        multiRunResultHistograms[(i*4)+2]->Draw();
        canvasArray[i]->cd(4);
        multiRunResultHistograms[(i*4)+3]->Draw();
    }
}

//fills the histograms, first resets histograms then fills the array with the lambda value fitted at run j
TH1D** Run::fillRunResultHistos(TH1D** multiRunResultHistograms)
{
    for(int i = 0; i < numElements; i++)
    {
        multiRunResultHistograms[(i*2)]->Reset("ICES");
        multiRunResultHistograms[(i*2)+1]->Reset("ICES");
        for(int j = 0; j < runs; j++)
        {
            multiRunResultHistograms[(i*2)]->Fill(batemanFitValues->GetAnHalfLife(j, i));
            multiRunResultHistograms[(i*2)+1]->Fill(integralFitValues->GetAnHalfLife(j, i));
        }
    }

    return multiRunResultHistograms;
}

//newB fills the single element histogram
TH1D** Run::fillRunResultHistosSingleElement(TH1D** multiRunResultHistogramsSingleElement)
{
    for(int i = 0; i < numElements; i++)
    {
        multiRunResultHistogramsSingleElement[(i*2)]->Reset("ICES");
        multiRunResultHistogramsSingleElement[(i*2)+1]->Reset("ICES");
        for(int j = 0; j < runs; j++)
        {
            multiRunResultHistogramsSingleElement[(i*2)]->Fill(singleBatemanFitValues->GetAnHalfLife(j, i, i));
            multiRunResultHistogramsSingleElement[(i*2)+1]->Fill(singleIntegralFitValues->GetAnHalfLife(j, i, i));
        }
    }

    return multiRunResultHistogramsSingleElement;
}

//fill graphs for the result of the multiple runs, used for events changing between runs (dynamic)
void Run::genGraphsEventChange()
{
    for(int j = 0; j < numElements; j++)
    {
        totalBatemanGraphs[j]= new TGraphErrors(runs, eventsXAxis, batemanFitValues->GetHalfLifeArr(j), zero, batemanFitValues->GetHalfLifeErrorArr(j));
        totalBatemanGraphs[j]->GetXaxis()->SetTitle("Events");
        totalBatemanGraphs[j]->GetYaxis()->SetTitle("Bateman Fit(S)");
        totalBatemanGraphs[j]->SetTitle((elementNameStrs[j] + " Bateman Fit(S)").c_str());
        totalBatemanGraphs[j]->SetName((elementNameStrs[j] + " Bateman_Fit(S)").c_str());

        totalIntegralGraphs[j] = new TGraphErrors(runs, eventsXAxis, integralFitValues->GetHalfLifeArr(j), zero, integralFitValues->GetHalfLifeErrorArr(j));
        totalIntegralGraphs[j]->GetXaxis()->SetTitle("Events");
        totalIntegralGraphs[j]->GetYaxis()->SetTitle("Integral Fit(S)");
        totalIntegralGraphs[j]->SetTitle((elementNameStrs[j] + " Integral Fit(S)").c_str());
        totalIntegralGraphs[j]->SetName((elementNameStrs[j] + " Integral_Fit(S)").c_str());
    }
}

//fill graphs for the result of the multiple runs, used for no change between runs (dynamic)
void Run::genGraphsNoChange() 
{
    for(int j = 0; j < numElements; j++)
    {
        totalBatemanGraphs[j]= new TGraphErrors(runs, runsXAxis, batemanFitValues->GetHalfLifeArr(j), zero, batemanFitValues->GetHalfLifeErrorArr(j));
        totalBatemanGraphs[j]->GetXaxis()->SetTitle("Runs");
        totalBatemanGraphs[j]->GetYaxis()->SetTitle("Bateman Fit(S)");
        totalBatemanGraphs[j]->SetTitle((elementNameStrs[j] + " bateman Fit(S)").c_str());
        totalBatemanGraphs[j]->SetName((elementNameStrs[j] + " bateman_Fit(S)").c_str());

        totalIntegralGraphs[j] = new TGraphErrors(runs, runsXAxis, integralFitValues->GetHalfLifeArr(j), zero, integralFitValues->GetHalfLifeErrorArr(j));
        totalIntegralGraphs[j]->GetXaxis()->SetTitle("Runs");
        totalIntegralGraphs[j]->GetYaxis()->SetTitle("Integral Fit(S)");
        totalIntegralGraphs[j]->SetTitle((elementNameStrs[j] + " Integral Fit(S)").c_str());
        totalIntegralGraphs[j]->SetName((elementNameStrs[j] + " Integral_Fit(S)").c_str());
    }
}

//newB generates the graphs for the single fit function data (dynamic)
void Run::genGraphsNoChangeSingleElement()
{
    for(int j = 0; j < numElements; j++)
    {
        singleBatemanGraphs[j]= new TGraphErrors(runs, runsXAxis, singleBatemanFitValues->GetHalfLifeArr(j, j), zero, singleBatemanFitValues->GetHalfLifeErrorArr(j, j));
        singleBatemanGraphs[j]->GetXaxis()->SetTitle("Runs");
        singleBatemanGraphs[j]->GetYaxis()->SetTitle("Bateman Fit(S)");
        singleBatemanGraphs[j]->SetTitle((elementNameStrs[j] + " bateman Fit Single Element(S)").c_str());
        singleBatemanGraphs[j]->SetName((elementNameStrs[j] + " bateman_Fit_Single_Element(S)").c_str());

        singleIntegralGraphs[j] = new TGraphErrors(runs, runsXAxis, singleIntegralFitValues->GetHalfLifeArr(j, j), zero, singleIntegralFitValues->GetHalfLifeErrorArr(j, j));
        singleIntegralGraphs[j]->GetXaxis()->SetTitle("Runs");
        singleIntegralGraphs[j]->GetYaxis()->SetTitle("Integral Fit Single Element(S)");
        singleIntegralGraphs[j]->SetTitle((elementNameStrs[j] + " Integral Fit Single Element(S)").c_str());
        singleIntegralGraphs[j]->SetName((elementNameStrs[j] + " Integral_Fit_Single_Element(S)").c_str());
    }
}

//returns maximum value in array
Double_t Run::getMaxElement(Double_t* arr)
{
    Double_t hold = arr[0];
    for(int i = 1; i < runs; i++)
    {
        if(hold < arr[i])
        {
            hold  = arr[i];
        }
    }
    return hold;
}

//returns minimum value in array
Double_t Run::getMinElement(Double_t* arr)
{
    Double_t hold = arr[0];
    for(int i = 1; i < runs; i++)
    {
        if(hold > arr[i])
        {
            hold  = arr[i];
        }
    }
    return hold;
}

//does run for no change and puts data in respective array
void Run::runNoChange(Int_t cycleIndex)
{
    //dynamic array
    ChainFitValues* tempFitParameters;
    SingleElementFitValues* singleTempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    for(int j = 0; j < runs; j++)
    {
        //generates random data and fits it. Then extract the fit parametes
        element->fitHistos(cycleIndex, j);
        //gets total bateman fit parameters
        tempFitParameters = element->getBatemanFitParameters();
        
        //move total bateman parameters into respective class
        for(int i = 0; i < numElements; i++)
        {
            tempN0 = tempFitParameters->GetAnN0(i);
            tempN0Error = tempFitParameters->GetAnN0Error(i);
            tempHalfLife = tempFitParameters->GetAnHalfLife(i);
            tempHalfLifeError = tempFitParameters->GetAnHalfLifeError(i);

            batemanFitValues->SetAnN0(j, i, tempN0);
            batemanFitValues->SetAnN0Error(j, i, tempN0Error);
            batemanFitValues->SetAnHalfLife(j, i, tempHalfLife);
            batemanFitValues->SetAnHalfLifeError(j, i, tempHalfLifeError);
        }

        //gets total integral fit parameters
        tempFitParameters = element->getIntegralFitParameters();

        //move total integral parameters into respective class
        for(int i = 0; i < numElements; i++)
        {
            tempN0 = tempFitParameters->GetAnN0(i);
            tempN0Error = tempFitParameters->GetAnN0Error(i);
            tempHalfLife = tempFitParameters->GetAnHalfLife(i);
            tempHalfLifeError = tempFitParameters->GetAnHalfLifeError(i);

            integralFitValues->SetAnN0(j, i, tempN0);
            integralFitValues->SetAnN0Error(j, i, tempN0Error);
            integralFitValues->SetAnHalfLife(j, i, tempHalfLife);
            integralFitValues->SetAnHalfLifeError(j, i, tempHalfLifeError);
        }

        //get single bateman fit parameters
        singleTempFitParameters = element->getSingleBatemanFitParameters();

        //move single bateman fit parameters into respective class
        for(int i = 0; i < numElements; i++)
        {
            for(int subIndex = 0; subIndex < i+1; subIndex++)
            {
                tempN0 = singleTempFitParameters->GetAnN0(i, subIndex);
                singleBatemanFitValues->SetAnN0(j, i, subIndex, tempN0);
                tempN0Error = singleTempFitParameters->GetAnN0Error(i, subIndex);
                singleBatemanFitValues->SetAnN0Error(j, i, subIndex, tempN0Error);
                tempHalfLife = singleTempFitParameters->GetAnHalfLife(i, subIndex);
                singleBatemanFitValues->SetAnHalfLife(j, i, subIndex, tempHalfLife);
                tempHalfLifeError = singleTempFitParameters->GetAnHalfLifeError(i, subIndex);
                singleBatemanFitValues->SetAnHalfLifeError(j, i, subIndex, tempHalfLifeError);
            }
        }

        //get single integral fit parameters
        singleTempFitParameters = element->getSingleIntegralFitParameters();

        //move single integral fit parameters into respective class
        for(int i = 0; i < numElements; i++)
        {
            for(int subIndex = 0; subIndex < i+1; subIndex++)
            {
                tempN0 = singleTempFitParameters->GetAnN0(i, subIndex);
                singleIntegralFitValues->SetAnN0(j, i, subIndex, tempN0);
                tempN0Error = singleTempFitParameters->GetAnN0Error(i, subIndex);
                singleIntegralFitValues->SetAnN0Error(j, i, subIndex, tempN0Error);
                tempHalfLife = singleTempFitParameters->GetAnHalfLife(i, subIndex);
                singleIntegralFitValues->SetAnHalfLife(j, i, subIndex, tempHalfLife);
                tempHalfLifeError = singleTempFitParameters->GetAnHalfLifeError(i, subIndex);
                singleIntegralFitValues->SetAnHalfLifeError(j, i, subIndex, tempHalfLifeError);
            }
        }
    }
}


//fits the singular histogram and puts the data of the singlular fit into the arrays (dynamic)
void Run::runNoChangeGenOnce(Int_t cycleIndex, Int_t runIndex)
{
    ChainFitValues* tempFitParameters;
    SingleElementFitValues* singleTempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    element->fitDataGenOnce(cycleIndex, runIndex);

    tempFitParameters = element->getBatemanFitParameters();
        
        //move total bateman parameters into respective class
        for(int i = 0; i < numElements; i++)
        {
            tempN0 = tempFitParameters->GetAnN0(i);
            tempN0Error = tempFitParameters->GetAnN0Error(i);
            tempHalfLife = tempFitParameters->GetAnHalfLife(i);
            tempHalfLifeError = tempFitParameters->GetAnHalfLifeError(i);

            batemanFitValues->SetAnN0(0, i, tempN0);
            batemanFitValues->SetAnN0Error(0, i, tempN0Error);
            batemanFitValues->SetAnHalfLife(0, i, tempHalfLife);
            batemanFitValues->SetAnHalfLifeError(0, i, tempHalfLifeError);
        }

        //gets total integral fit parameters
        tempFitParameters = element->getIntegralFitParameters();

        //move total integral parameters into respective class
        for(int i = 0; i < numElements; i++)
        {
            tempN0 = tempFitParameters->GetAnN0(i);
            tempN0Error = tempFitParameters->GetAnN0Error(i);
            tempHalfLife = tempFitParameters->GetAnHalfLife(i);
            tempHalfLifeError = tempFitParameters->GetAnHalfLifeError(i);

            integralFitValues->SetAnN0(0, i, tempN0);
            integralFitValues->SetAnN0Error(0, i, tempN0Error);
            integralFitValues->SetAnHalfLife(0, i, tempHalfLife);
            integralFitValues->SetAnHalfLifeError(0, i, tempHalfLifeError);
        }

        //get single bateman fit parameters
        singleTempFitParameters = element->getSingleBatemanFitParameters();

        //move single bateman fit parameters into respective class
        for(int i = 0; i < numElements; i++)
        {
            for(int subIndex = 0; subIndex < i+1; subIndex++)
            {
                tempN0 = singleTempFitParameters->GetAnN0(i, subIndex);
                tempN0Error = singleTempFitParameters->GetAnN0Error(i, subIndex);
                tempHalfLife = singleTempFitParameters->GetAnHalfLife(i, subIndex);
                tempHalfLifeError = singleTempFitParameters->GetAnHalfLifeError(i, subIndex);

                singleBatemanFitValues->SetAnN0(0, i, subIndex, tempN0);
                singleBatemanFitValues->SetAnN0Error(0, i, subIndex, tempN0Error);
                singleBatemanFitValues->SetAnHalfLife(0, i, subIndex, tempHalfLife);
                singleBatemanFitValues->SetAnHalfLifeError(0, i, subIndex, tempHalfLifeError);
            }
        }

        //get single integral fit parameters
        singleTempFitParameters = element->getSingleIntegralFitParameters();

        //move single integral fit parameters into respective class
        for(int i = 0; i < numElements; i++)
        {
            for(int subIndex = 0; subIndex < i+1; subIndex++)
            {
                tempN0 = singleTempFitParameters->GetAnN0(i, subIndex);
                singleIntegralFitValues->SetAnN0(0, i, subIndex, tempN0);
                tempN0Error = singleTempFitParameters->GetAnN0Error(i, subIndex);
                singleIntegralFitValues->SetAnN0Error(0, i, subIndex, tempN0Error);
                tempHalfLife = singleTempFitParameters->GetAnHalfLife(i, subIndex);
                singleIntegralFitValues->SetAnHalfLife(0, i, subIndex, tempHalfLife);
                tempHalfLifeError = singleTempFitParameters->GetAnHalfLifeError(i, subIndex);
                singleIntegralFitValues->SetAnHalfLifeError(0, i, subIndex, tempHalfLifeError);
            }
        }
}

#endif