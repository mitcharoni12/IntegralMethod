#ifndef RUN_H
#define RUN_H

#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cstdio>
#include <sys/stat.h>

#include "ElementFit.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "FitParameterStore.h"

using namespace std;

class Run{
    private:
        Int_t runs, numElements, eventDecrement;
        ElementFit* element; 
        Double_t** regularFitErrors, **integralFitErrors, **regularFitValue, **integralFitValue, **initRegularValue, **initIntegralValue,
                 **regularFitErrorsSingleElement, **integralFitErrorsSingleElement, **regularFitValueSingleElement, **integralFitValueSingleElement, **initRegularValueSingleElement, **initIntegralValueSingleElement;
        Double_t* eventsXAxis, *runsXAxis;
        string* elementNameStrs;
        TH3D* correlationHisto;
        TH1D** multiRunResultHistograms;
        TCanvas* testCan;
        TCanvas** multiRunResultCanvas;
        //helper functions
        void createCorrelationHistoEventChange();
        void createCorrelationHistoNoChange();
        Double_t getMaxElement(Double_t* arr);
        Double_t getMinElement(Double_t* arr);
        void runEventChange();
        void genIntegralMeanGraphs();
        void genRegularMeanGraphs();
    public:
        Run(Int_t runs, Int_t eventDecrement, ElementFit* element, string* elementNameStrs);
        ~Run();
        Double_t** runNoChangeGenOnce();
        TGraph** genGraphsEventChange();
        TGraph** genGraphsNoChange();
        TGraph** genGraphsNoChangeSingleElement();
        TH1D** createRunResultHistos();
        TH1D** createRunResultHistosSingleElements();
        TH1D** fillRunResultHistos(TH1D** multiRunResultHistograms);
        TH1D** fillRunResultHistosSingleElement(TH1D** multiRunResultHistogramsSingleElement);
        void displayMultiRunResultGraphs(TCanvas** canvasArray, TGraph** multiRunResultGraph);
        void displayCorrelationHisto(TCanvas* canvas);
        void displayMultiRunResultHistos(TCanvas** canvasArray, TH1D** multiRunResultHistograms);
        void runNoChange();
        //getter function
        TH1D** getMultiRunResultHistos(){return multiRunResultHistograms;}
        string* elementStringNames(){return elementNameStrs;}
        //setter function
        void setNumRuns(Int_t numRuns){this->runs = numRuns;}
};

Run::Run(Int_t runs, Int_t eventDecrement, ElementFit* element, string* elementNameStrs)
{
    //variable parameter setting
    this->runs = runs;
    this->eventDecrement = eventDecrement;
    this->element = element;
    this->elementNameStrs = elementNameStrs;
    numElements = element->getNumElements();

    //dynamically allocating required arrays and root variables
    //dynamic array
    regularFitErrors = new Double_t* [numElements];
    integralFitErrors = new Double_t* [numElements];
    regularFitValue = new Double_t* [numElements];
    integralFitValue = new Double_t* [numElements];
    initRegularValue = new Double_t* [numElements];
    initIntegralValue = new Double_t* [numElements];
    regularFitErrorsSingleElement = new Double_t* [numElements];
    integralFitErrorsSingleElement = new Double_t* [numElements];
    regularFitValueSingleElement = new Double_t* [numElements];
    integralFitValueSingleElement = new Double_t* [numElements];
    initRegularValueSingleElement = new Double_t* [numElements];
    initIntegralValueSingleElement = new Double_t* [numElements];
    eventsXAxis = new Double_t [runs];
    runsXAxis = new Double_t [runs];
    //testCan = new TCanvas("testCan", "testCan", 500, 500);
    for(int i = 0; i < numElements; i++)
    {
        regularFitErrors[i] = new Double_t [runs];
        integralFitErrors[i] = new Double_t [runs];
        regularFitValue[i] = new Double_t [runs];
        integralFitValue[i] = new Double_t [runs];
        initRegularValue[i] = new Double_t [runs];
        initIntegralValue[i] = new Double_t [runs];
        regularFitErrorsSingleElement[i] = new Double_t [runs];
        integralFitErrorsSingleElement[i] = new Double_t [runs];
        regularFitValueSingleElement[i] = new Double_t [runs];
        integralFitValueSingleElement[i] = new Double_t [runs];
        initRegularValueSingleElement[i] = new Double_t [runs];
        initIntegralValueSingleElement[i] = new Double_t [runs];
    }
    for(int i = 0; i < runs; i++)
    {
        runsXAxis[i] = i+1;
    }
}

Run::~Run()
{
    for(int i = 0; i < element->getNumElements(); i++)
    {
        delete [] regularFitErrors[i];
        delete [] integralFitErrors[i];
        delete [] regularFitValue[i];
        delete [] integralFitValue[i];
        delete [] initRegularValue[i];
        delete [] initIntegralValue[i];
        delete [] regularFitErrorsSingleElement[i];
        delete [] integralFitErrorsSingleElement[i];
        delete [] regularFitValueSingleElement[i];
        delete [] integralFitValueSingleElement[i];
        delete [] initRegularValueSingleElement[i];
        delete [] initIntegralValueSingleElement[i];
    }
    delete [] regularFitErrors;
    delete [] integralFitErrors;
    delete [] regularFitValue;
    delete [] integralFitValue;
    delete [] initRegularValue;
    delete [] initIntegralValue;
    delete [] regularFitErrorsSingleElement;
    delete [] integralFitErrorsSingleElement;
    delete [] regularFitValueSingleElement;
    delete [] integralFitValueSingleElement;
    delete [] initRegularValueSingleElement;
    delete [] initIntegralValueSingleElement;
    delete [] eventsXAxis;
    delete [] runsXAxis;
}

//(WIP) creates and fills the correlation canvas for runs with event changes
void Run::createCorrelationHistoEventChange()
{
    correlationHisto = new TH3D("CorrelationHistoEventChange", "CorrelationHistoEventChange", 100, getMinElement(regularFitErrors[0]), getMaxElement(regularFitErrors[0]), 100, getMinElement(integralFitErrors[0]), getMaxElement(integralFitErrors[0]), 100, getMinElement(eventsXAxis), getMaxElement(eventsXAxis));
    correlationHisto->GetXaxis()->SetTitle("Reg Errors");
    correlationHisto->GetXaxis()->SetNdivisions(5);
    correlationHisto->GetXaxis()->CenterTitle(kTRUE);
    correlationHisto->GetXaxis()->SetTitleColor(kRed);
    correlationHisto->GetYaxis()->SetTitle("Integral Errors");
    correlationHisto->GetYaxis()->SetNdivisions(5);
    correlationHisto->GetYaxis()->CenterTitle(kTRUE);
    correlationHisto->GetYaxis()->SetTitleColor(kRed);
    correlationHisto->GetZaxis()->SetTitle("Events");
    correlationHisto->GetZaxis()->SetNdivisions(5);
    correlationHisto->GetZaxis()->CenterTitle(kTRUE);
    correlationHisto->GetZaxis()->SetTitleColor(kRed);
    for(int i = 0; i < runs; i++)
    {
        correlationHisto->Fill(regularFitErrors[0][i], integralFitErrors[0][i], eventsXAxis[i]);
    }
}

//(WIP) creates and fills the correlation canvas for runs with no changes
void Run::createCorrelationHistoNoChange()
{
    correlationHisto = new TH3D("CorrelationHistoNoChange", "CorrelationHistoNoChange", 100, getMinElement(regularFitErrors[0]), getMaxElement(regularFitErrors[0]), 100, getMinElement(integralFitErrors[0]), getMaxElement(integralFitErrors[0]), 100, 0, runs);
    correlationHisto->GetXaxis()->SetTitle("Reg Errors");
    correlationHisto->GetXaxis()->SetNdivisions(5);
    correlationHisto->GetXaxis()->CenterTitle(kTRUE);
    correlationHisto->GetXaxis()->SetTitleColor(kRed);
    correlationHisto->GetYaxis()->SetTitle("Integral Errors");
    correlationHisto->GetYaxis()->SetNdivisions(5);
    correlationHisto->GetYaxis()->CenterTitle(kTRUE);
    correlationHisto->GetYaxis()->SetTitleColor(kRed);
    correlationHisto->GetZaxis()->SetTitle("runs");
    correlationHisto->GetZaxis()->SetNdivisions(5);
    correlationHisto->GetZaxis()->CenterTitle(kTRUE);
    correlationHisto->GetZaxis()->SetTitleColor(kRed);
    for(int i = 0; i < runs; i++)
    {
        correlationHisto->Fill(regularFitErrors[0][i], integralFitErrors[0][i], runsXAxis[i]);
    }
}

//creates histograms for the run result of the total functions (dynamic)
TH1D** Run::createRunResultHistos()
{
    TH1D** multiRunResultHistograms = new TH1D* [numElements*2];
    Double_t parameterValue;

    //creates the histograms
    for(int i = 0; i < numElements; i++)
    {
        parameterValue = TMath::LogE()/(element->getElementParameters(i));
        multiRunResultHistograms[(i*2)] = new TH1D((elementNameStrs[i] + " fitResultRegHisto").c_str(), (elementNameStrs[i] + " fitResultRegHisto").c_str(), 500, parameterValue*0, parameterValue*2.5);
        multiRunResultHistograms[(i*2)+1] = new TH1D((elementNameStrs[i] + " fitResultIntegralHisto").c_str(), (elementNameStrs[i] + " fitResultIntegralHisto").c_str(), 500, parameterValue*0, parameterValue*2.5);
    }

    return multiRunResultHistograms;
}

//newB creates histograms for the run result of the single element functions (dynamcic)
TH1D** Run::createRunResultHistosSingleElements()
{
    TH1D** multiRunResultHistosSingleElement = new TH1D* [numElements*2];
    Double_t parameterValue;

    //creates the histograms
    for(int i = 0; i < numElements; i++)
    {
        parameterValue = TMath::LogE()/(element->getElementParameters(i));
        multiRunResultHistosSingleElement[(i*2)] = new TH1D((elementNameStrs[i] + " fitResultRegHistoSingleElement").c_str(), (elementNameStrs[i] + " fitResultRegHistoSingleElement").c_str(), 500, parameterValue*0, parameterValue*2.5);
        multiRunResultHistosSingleElement[(i*2)+1] = new TH1D((elementNameStrs[i] + " fitResultIntegralHistoSingleElement").c_str(), (elementNameStrs[i] + " fitResultIntegralHistoSingleElement").c_str(), 500, parameterValue*0, parameterValue*2.5);
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
void Run::displayMultiRunResultGraphs(TCanvas** canvasArray, TGraph** multiRunResultGraph)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(1);
        multiRunResultGraph[(i*4)]->Draw();
        canvasArray[i]->cd(2);
        multiRunResultGraph[(i*4)+1]->Draw();
        canvasArray[i]->cd(3);
        multiRunResultGraph[(i*4)+2]->Draw();
        canvasArray[i]->cd(4);
        multiRunResultGraph[(i*4)+3]->Draw();
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
            multiRunResultHistograms[(i*2)]->Fill(regularFitValue[i][j]);
            multiRunResultHistograms[(i*2)+1]->Fill(integralFitValue[i][j]);
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
            multiRunResultHistogramsSingleElement[(i*2)]->Fill(regularFitValueSingleElement[i][j]);
            multiRunResultHistogramsSingleElement[(i*2)+1]->Fill(integralFitValueSingleElement[i][j]);
        }
    }

    return multiRunResultHistogramsSingleElement;
}

//fill graphs for the result of the multiple runs, used for events changing between runs (dynamic)
TGraph** Run::genGraphsEventChange()
{
    TGraph** multiRunResultGraph = new TGraph* [numElements*6];

    for(int j = 0; j < numElements; j++)
    {
        multiRunResultGraph[j*6] = new TGraph(runs, eventsXAxis, regularFitErrors[j]);
        multiRunResultGraph[j*6]->GetXaxis()->SetTitle(("events" + to_string((j*6)+1)).c_str());
        multiRunResultGraph[j*6]->GetYaxis()->SetTitle(("Regular Fit Error" + to_string((j*6)+1)).c_str());
        multiRunResultGraph[j*6]->SetTitle((elementNameStrs[j] + "Regular Fit Error").c_str());
        multiRunResultGraph[j*6]->SetName((elementNameStrs[j] + "Regular_Fit_Error").c_str());

        multiRunResultGraph[(j*6)+1] = new TGraph(runs, eventsXAxis, integralFitErrors[j]);
        multiRunResultGraph[(j*6)+1]->GetXaxis()->SetTitle(("events" + to_string((j*6)+2)).c_str());
        multiRunResultGraph[(j*6)+1]->GetYaxis()->SetTitle(("Integral Fit Error" + to_string((j*6)+2)).c_str());
        multiRunResultGraph[(j*6)+1]->SetTitle((elementNameStrs[j] + "Integral Fit Error").c_str());
        multiRunResultGraph[(j*6)+1]->SetName((elementNameStrs[j] + "Integral_Fit_Error").c_str());

        multiRunResultGraph[(j*6)+2] = new TGraph(runs, eventsXAxis, regularFitValue[j]);
        multiRunResultGraph[(j*6)+2]->GetXaxis()->SetTitle(("events" + to_string((j*6)+3)).c_str());
        multiRunResultGraph[(j*6)+2]->GetYaxis()->SetTitle(("Regular Fit Value" + to_string((j*6)+3)).c_str());
        multiRunResultGraph[(j*6)+2]->SetTitle((elementNameStrs[j] + "Regular Fit Value").c_str());
        multiRunResultGraph[(j*6)+2]->SetName((elementNameStrs[j] + "Regular_Fit_Value").c_str());

        multiRunResultGraph[(j*6)+3] = new TGraph(runs, eventsXAxis, integralFitValue[j]);
        multiRunResultGraph[(j*6)+3]->GetXaxis()->SetTitle(("events" + to_string((j*6)+4)).c_str());
        multiRunResultGraph[(j*6)+3]->GetYaxis()->SetTitle(("Integral Fit Value" + to_string((j*6)+4)).c_str());
        multiRunResultGraph[(j*6)+3]->SetTitle((elementNameStrs[j] + "Integral Fit Value").c_str());
        multiRunResultGraph[(j*6)+3]->SetName((elementNameStrs[j] + "Integral_Fit_Value").c_str());

        multiRunResultGraph[(j*6)+4] = new TGraph(runs, eventsXAxis, initRegularValue[j]);
        multiRunResultGraph[(j*6)+4]->GetXaxis()->SetTitle(("events" + to_string((j*6)+5)).c_str());
        multiRunResultGraph[(j*6)+4]->GetYaxis()->SetTitle(("Init Regular Value" + to_string((j*6)+5)).c_str());
        multiRunResultGraph[(j*6)+4]->SetTitle((elementNameStrs[j] + "Init Regular value").c_str());
        multiRunResultGraph[(j*6)+4]->SetName((elementNameStrs[j] + "Init_Regular_value").c_str());

        multiRunResultGraph[(j*6)+5] = new TGraph(runs, eventsXAxis, initIntegralValue[j]);
        multiRunResultGraph[(j*6)+5]->GetXaxis()->SetTitle(("events" + to_string((j*6)+6)).c_str());
        multiRunResultGraph[(j*6)+5]->GetYaxis()->SetTitle(("Init Integral Value" + to_string((j*6)+6)).c_str());
        multiRunResultGraph[(j*6)+5]->SetTitle((elementNameStrs[j] + "Init Integral value").c_str());
        multiRunResultGraph[(j*6)+5]->SetName((elementNameStrs[j] + "Init_Integral_value").c_str());
    }
    return multiRunResultGraph;
}

//fill graphs for the result of the multiple runs, used for no change between runs (dynamic)
TGraph** Run::genGraphsNoChange() 
{
    TGraph** multiRunResultGraph = new TGraph* [numElements*6];

    for(int j = 0; j < numElements; j++)
    {
        multiRunResultGraph[j*6] = new TGraph(runs, runsXAxis, regularFitErrors[j]);
        multiRunResultGraph[j*6]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+1)).c_str());
        multiRunResultGraph[j*6]->GetYaxis()->SetTitle(("Regular Fit Error" + to_string((j*6)+1)).c_str());
        multiRunResultGraph[j*6]->SetTitle((elementNameStrs[j] + "Regular Fit Error").c_str());
        multiRunResultGraph[j*6]->SetName((elementNameStrs[j] + "Regular_Fit_Error").c_str());

        multiRunResultGraph[(j*6)+1] = new TGraph(runs, runsXAxis, integralFitErrors[j]);
        multiRunResultGraph[(j*6)+1]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+2)).c_str());
        multiRunResultGraph[(j*6)+1]->GetYaxis()->SetTitle(("Integral Fit Error" + to_string((j*6)+2)).c_str());
        multiRunResultGraph[(j*6)+1]->SetTitle((elementNameStrs[j] + "Integral Fit Error").c_str());
        multiRunResultGraph[(j*6)+1]->SetName((elementNameStrs[j] + "Integral_Fit_Error").c_str());

        multiRunResultGraph[(j*6)+2] = new TGraph(runs, runsXAxis, regularFitValue[j]);
        multiRunResultGraph[(j*6)+2]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+3)).c_str());
        multiRunResultGraph[(j*6)+2]->GetYaxis()->SetTitle(("Regular Fit Value" + to_string((j*6)+3)).c_str());
        multiRunResultGraph[(j*6)+2]->SetTitle((elementNameStrs[j] + "Regular Fit Value").c_str());
        multiRunResultGraph[(j*6)+2]->SetName((elementNameStrs[j] + "Regular_Fit_Value").c_str());

        multiRunResultGraph[(j*6)+3] = new TGraph(runs, runsXAxis, integralFitValue[j]);
        multiRunResultGraph[(j*6)+3]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+4)).c_str());
        multiRunResultGraph[(j*6)+3]->GetYaxis()->SetTitle(("Integral Fit Value" + to_string((j*6)+4)).c_str());
        multiRunResultGraph[(j*6)+3]->SetTitle((elementNameStrs[j] + "Integral Fit Value").c_str());
        multiRunResultGraph[(j*6)+3]->SetName((elementNameStrs[j] + "Integral_Fit_Value").c_str());

        multiRunResultGraph[(j*6)+4] = new TGraph(runs, runsXAxis, initRegularValue[j]);
        multiRunResultGraph[(j*6)+4]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+5)).c_str());
        multiRunResultGraph[(j*6)+4]->GetYaxis()->SetTitle(("Init Regular Value" + to_string((j*6)+5)).c_str());
        multiRunResultGraph[(j*6)+4]->SetTitle((elementNameStrs[j] + "Init Regular value").c_str());
        multiRunResultGraph[(j*6)+4]->SetName((elementNameStrs[j] + "Init_Regular_value").c_str());

        multiRunResultGraph[(j*6)+5] = new TGraph(runs, runsXAxis, initIntegralValue[j]);
        multiRunResultGraph[(j*6)+5]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+6)).c_str());
        multiRunResultGraph[(j*6)+5]->GetYaxis()->SetTitle(("Init Integral Value" + to_string((j*6)+6)).c_str());
        multiRunResultGraph[(j*6)+5]->SetTitle((elementNameStrs[j] + "Init Integral value").c_str());
        multiRunResultGraph[(j*6)+5]->SetName((elementNameStrs[j] + "Init_Integral_value").c_str());
    }
    return multiRunResultGraph;
}

//newB generates the graphs for the single fit function data (dynamic)
TGraph** Run::genGraphsNoChangeSingleElement()
{
    TGraph** multiRunResultGraphsSingleElement = new TGraph* [numElements*6];

    for(int j = 0; j < numElements; j++)
    {
        multiRunResultGraphsSingleElement[j*6] = new TGraph(runs, runsXAxis, regularFitErrorsSingleElement[j]);
        multiRunResultGraphsSingleElement[j*6]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+1)).c_str());
        multiRunResultGraphsSingleElement[j*6]->GetYaxis()->SetTitle(("Regular Fit Error" + to_string((j*6)+1)).c_str());
        multiRunResultGraphsSingleElement[j*6]->SetTitle((elementNameStrs[j] + "Regular Fit Error Single Element").c_str());
        multiRunResultGraphsSingleElement[j*6]->SetName((elementNameStrs[j] + "Regular_Fit_Error_Single_Element").c_str());

        multiRunResultGraphsSingleElement[(j*6)+1] = new TGraph(runs, runsXAxis, integralFitErrorsSingleElement[j]);
        multiRunResultGraphsSingleElement[(j*6)+1]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+2)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+1]->GetYaxis()->SetTitle(("Integral Fit Error" + to_string((j*6)+2)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+1]->SetTitle((elementNameStrs[j] + "Integral Fit Error Single Element").c_str());
        multiRunResultGraphsSingleElement[(j*6)+1]->SetName((elementNameStrs[j] + "Integral_Fit_Error_Single_Element").c_str());

        multiRunResultGraphsSingleElement[(j*6)+2] = new TGraph(runs, runsXAxis, regularFitValueSingleElement[j]);
        multiRunResultGraphsSingleElement[(j*6)+2]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+3)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+2]->GetYaxis()->SetTitle(("Regular Fit Value" + to_string((j*6)+3)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+2]->SetTitle((elementNameStrs[j] + "Regular Fit Value Single Element").c_str());
        multiRunResultGraphsSingleElement[(j*6)+2]->SetName((elementNameStrs[j] + "Regular_Fit_Value_Single_Element").c_str());

        multiRunResultGraphsSingleElement[(j*6)+3] = new TGraph(runs, runsXAxis, integralFitValueSingleElement[j]);
        multiRunResultGraphsSingleElement[(j*6)+3]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+4)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+3]->GetYaxis()->SetTitle(("Integral Fit Value" + to_string((j*6)+4)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+3]->SetTitle((elementNameStrs[j] + "Integral Fit Value Single Element").c_str());
        multiRunResultGraphsSingleElement[(j*6)+3]->SetName((elementNameStrs[j] + "Integral_Fit_Value_Single_Element").c_str());

        multiRunResultGraphsSingleElement[(j*6)+4] = new TGraph(runs, runsXAxis, initRegularValueSingleElement[j]);
        multiRunResultGraphsSingleElement[(j*6)+4]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+5)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+4]->GetYaxis()->SetTitle(("Init Regular Value" + to_string((j*6)+5)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+4]->SetTitle((elementNameStrs[j] + "Init Regular value Single Element").c_str());
        multiRunResultGraphsSingleElement[(j*6)+4]->SetName((elementNameStrs[j] + "Init_Regular_value_Single_Element").c_str());

        multiRunResultGraphsSingleElement[(j*6)+5] = new TGraph(runs, runsXAxis, initIntegralValueSingleElement[j]);
        multiRunResultGraphsSingleElement[(j*6)+5]->GetXaxis()->SetTitle(("Runs" + to_string((j*6)+6)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+5]->GetYaxis()->SetTitle(("Init Integral Value" + to_string((j*6)+6)).c_str());
        multiRunResultGraphsSingleElement[(j*6)+5]->SetTitle((elementNameStrs[j] + "Init Integral value Single Element").c_str());
        multiRunResultGraphsSingleElement[(j*6)+5]->SetName((elementNameStrs[j] + "Init_Integral_value_Single_Element").c_str());
    }
    return multiRunResultGraphsSingleElement;
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

//does run for changing events and puts data in respective array
void Run::runEventChange()
{
    FitParameterStore* fitParameters = element->getTotalFitParameters();
    for(int j = 0; j < runs; j++)
    {
        //generates random data and fits it. Then extract the fit parametes
        element->fitData();
        //takes fit parameters and puts them in their respective arrays
        for(int i = 0; i < numElements; i++)
        {
            eventsXAxis[j] = ((Double_t) element->getNumEvents());
            regularFitErrors[i][j] = (fitParameters->getRegularHalfLifeError())[i];
            integralFitErrors[i][j] = (fitParameters->getIntegralHalfLifeError())[i];
            regularFitValue[i][j] = (fitParameters->getRegularHalfLife())[i];
            integralFitValue[i][j] = (fitParameters->getIntegralHalfLife())[i];
            initRegularValue[i][j] = (fitParameters->getRegularN0())[i];
            initIntegralValue[i][j] = (fitParameters->getIntegralN0())[i];
        }
        //changes number of events
        element->setNumEvents(element->getNumEvents()- eventDecrement);
    }
}

//does run for no change and puts data in respective array
void Run::runNoChange()
{
    //dynamic array
    FitParameterStore* fitParameters;
    FitParameterStore** singleElementFitParameters;
    for(int j = 0; j < runs; j++)
    {
        //generates random data and fits it. Then extract the fit parametes
        element->fitData();
        fitParameters = element->getTotalFitParameters();
        singleElementFitParameters = element->getSingleFitParameters();
        
        //takes fit parameters and puts them in their respective arrays
        for(int i = 0; i < numElements; i++)
        {
            regularFitErrors[i][j] = (fitParameters->getRegularHalfLifeError())[i];
            integralFitErrors[i][j] = (fitParameters->getIntegralHalfLifeError())[i];
            regularFitValue[i][j] = (fitParameters->getRegularHalfLife())[i];
            integralFitValue[i][j] = (fitParameters->getIntegralHalfLife())[i];
            initRegularValue[i][j] = (fitParameters->getRegularN0())[i];
            initIntegralValue[i][j] = (fitParameters->getIntegralN0())[i];
        }
        //newB gets fit parameters for the single element functions and puts them in their respective array
        for(int i = 0; i < numElements; i++)
        {
            regularFitErrorsSingleElement[i][j] = (singleElementFitParameters[i]->getRegularHalfLifeError())[i];
            integralFitErrorsSingleElement[i][j] = (singleElementFitParameters[i]->getIntegralHalfLifeError())[i];
            regularFitValueSingleElement[i][j] = (singleElementFitParameters[i]->getRegularHalfLife())[i];
            integralFitValueSingleElement[i][j] = (singleElementFitParameters[i]->getIntegralHalfLife())[i];
            initRegularValueSingleElement[i][j] = (singleElementFitParameters[i]->getRegularN0())[i];
            initIntegralValueSingleElement[i][j] = (singleElementFitParameters[i]->getIntegralN0())[i];
        }
    }
}

//fits the singular histogram and puts the data of the singlular fit into the arrays (dynamic)
Double_t** Run::runNoChangeGenOnce()
{
    FitParameterStore* fitParameters;
    FitParameterStore** singleFitParameters;
    Double_t** value = new Double_t* [12];

    Double_t* regularFitErrors = new Double_t [numElements];
    Double_t* integralFitErrors = new Double_t [numElements];
    Double_t* regularFitValue = new Double_t [numElements];
    Double_t* integralFitValue = new Double_t [numElements];
    Double_t* initRegularValue = new Double_t [numElements];
    Double_t* initIntegralValue = new Double_t [numElements];

    Double_t* regularFitErrorsSingleElement = new Double_t [numElements];
    Double_t* integralFitErrorsSingleElement = new Double_t [numElements];
    Double_t* regularFitValueSingleElement = new Double_t [numElements];
    Double_t* integralFitValueSingleElement= new Double_t [numElements];
    Double_t* initRegularValueSingleElement = new Double_t [numElements];
    Double_t* initIntegralValueSingleElement = new Double_t [numElements];

    element->fitDataGenOnce();
    fitParameters = element->getTotalFitParameters();
    singleFitParameters = element->getSingleFitParameters();
    //takes fit parameters and puts them in their respective arrays
    for(int i = 0; i < numElements; i++)
    {
        regularFitErrors[i] = (fitParameters->getRegularHalfLifeError())[i];
        integralFitErrors[i] = (fitParameters->getIntegralHalfLifeError())[i];
        regularFitValue[i] = (fitParameters->getRegularHalfLife())[i];
        integralFitValue[i] = (fitParameters->getIntegralHalfLife())[i];
        initRegularValue[i] = (fitParameters->getRegularN0())[i];
        initIntegralValue[i] = (fitParameters->getIntegralN0())[i];

        regularFitErrorsSingleElement[i] = (singleFitParameters[i]->getRegularHalfLifeError())[i];
        integralFitErrorsSingleElement[i] = (singleFitParameters[i]->getIntegralHalfLifeError())[i];
        regularFitValueSingleElement[i] = (singleFitParameters[i]->getRegularHalfLife())[i];
        integralFitValueSingleElement[i] = (singleFitParameters[i]->getIntegralHalfLife())[i];
        initRegularValueSingleElement[i] = (singleFitParameters[i]->getRegularN0())[i];
        initIntegralValueSingleElement[i] = (singleFitParameters[i]->getIntegralN0())[i];
    }
    
    value[0] = regularFitErrors;
    value[1] = integralFitErrors;
    value[2] = regularFitValue;
    value[3] = integralFitValue;
    value[4] = initRegularValue;
    value[5] = initIntegralValue;
    value[6] = regularFitErrorsSingleElement;
    value[7] = integralFitErrorsSingleElement;
    value[8] = regularFitValueSingleElement;
    value[9] = integralFitValueSingleElement;
    value[10] = initRegularValueSingleElement;
    value[11] = initIntegralValueSingleElement;

    return value;
}

#endif