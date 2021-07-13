#ifndef CYCLE_H
#define CYCLE_H

#include <iostream>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"

#include "Run.h"
#include "ElementFit.h"
#include "SingleChainRunFitValues.h"
#include "ChainRunFitValues.h"
#include "FitOption.h"

using namespace std;

class Cycle{
    private:
        FitOption* fitOptions;
        Int_t cycles, numElements, incChoice;
        Int_t* binNumArr;
        Double_t*timeArr;
        Run* decayChainRun;
        ElementFit* element;
        string* elementStrNames;
        ChainRunFitValues* batemanFitValues, *integralFitValues, *meanDifferenceValues;
        SingleChainRunFitValues* singleBatemanFitValues, *singleIntegralFitValues, *singleMeanDifferenceValues;
        TGraphErrors** batemanFitMeanGraphs, **integralFitMeanGraphs, **singleBatemanFitMeanGraphs
        ,**singleIntegralFitMeanGraphs, **fitDifferenceGraphs, **singleFitDifferenceGraphs;
        //TCanvas* test = new TCanvas("test", "test", 500, 500);
        //TCanvas* test2 = new TCanvas("test2", "test2", 500, 500);
        //helper functions
    public:
        Cycle(Run* decayChainRun, ElementFit* element);
        ~Cycle();
        void displayMeanDifferenceGraphs(TCanvas** canvasArr);
        void displayMeanSeperateGraphs(TCanvas** canvasArr);
        void genMeanDifferenceGraphsRebin();
        void genMeanDifferenceGraphsTimeDifference();
        void genSeperateMeanGraphsRebin();
        void genSeperateMeanGraphsTimeDifference();
        void genSingleMeanDifference();
        void runDifferenceMeanRebin();
        void runDifferenceMeanTimeDifference();
        void runSeperateMeanRebin();
        void runSeperateMeanTimeDifference();
        void runSeperateSingleGen();
};

Cycle::Cycle(Run* decayChainRun, ElementFit* element)
{
    this->fitOptions = element->getFitOptions();
    this->cycles = fitOptions->GetNumCycles();
    this->element = element;
    this->decayChainRun = decayChainRun;
    numElements = fitOptions->GetNumElements();
    binNumArr = fitOptions->GetBinNumArr();
    elementStrNames = fitOptions->GetElementNames();
    batemanFitValues = new ChainRunFitValues(numElements, cycles);
    integralFitValues = new ChainRunFitValues(numElements, cycles);
    meanDifferenceValues = new ChainRunFitValues(numElements, cycles);
    singleBatemanFitValues = new SingleChainRunFitValues(numElements, cycles);
    singleIntegralFitValues = new SingleChainRunFitValues(numElements, cycles);
    singleMeanDifferenceValues = new SingleChainRunFitValues(numElements, cycles);
    batemanFitMeanGraphs = new TGraphErrors* [numElements];
    integralFitMeanGraphs = new TGraphErrors* [numElements];
    singleBatemanFitMeanGraphs = new TGraphErrors* [numElements];
    singleIntegralFitMeanGraphs = new TGraphErrors* [numElements];
    fitDifferenceGraphs = new TGraphErrors* [numElements];
    singleFitDifferenceGraphs = new TGraphErrors* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        batemanFitMeanGraphs[i] = nullptr;
        integralFitMeanGraphs[i] = nullptr;
        singleBatemanFitMeanGraphs[i] = nullptr;
        singleIntegralFitMeanGraphs[i] = nullptr;
        fitDifferenceGraphs[i] = nullptr;
        singleFitDifferenceGraphs[i] = nullptr;
    }
    if(fitOptions->GetTimeShiftType() == 1)
    {
        timeArr = fitOptions->GetTimeFitEndArr();
    }else if(fitOptions->GetTimeShiftType() == 2)
    {
        timeArr = fitOptions->GetTimeFitStartArr();
    }
}

Cycle::~Cycle()
{
    delete [] batemanFitMeanGraphs;
    delete [] integralFitMeanGraphs;
    delete [] singleBatemanFitMeanGraphs;
    delete [] singleIntegralFitMeanGraphs;
    delete [] fitDifferenceGraphs;
    delete [] singleFitDifferenceGraphs;
    delete batemanFitValues;
    delete integralFitValues;
    delete meanDifferenceValues;
    delete singleBatemanFitValues;
    delete singleIntegralFitValues;
    delete singleMeanDifferenceValues;
}

//displays the graphs for the mean difference
void Cycle::displayMeanDifferenceGraphs(TCanvas** canvasArr)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(1);
        fitDifferenceGraphs[i]->Draw();
        canvasArr[i]->cd(2);
        singleFitDifferenceGraphs[i]->Draw();
    }
}

//displays the graphs for the seperate mean for time difference
void Cycle::displayMeanSeperateGraphs(TCanvas** canvasArr)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(1);
        batemanFitMeanGraphs[i]->Draw();
        canvasArr[i]->cd(2);
        integralFitMeanGraphs[i]->Draw();
        canvasArr[i]->cd(3);
        singleBatemanFitMeanGraphs[i]->Draw();
        canvasArr[i]->cd(4);
        singleIntegralFitMeanGraphs[i]->Draw();
    }
}

//creates and fills graphs for the difference in mean results for rebin
void Cycle::genMeanDifferenceGraphsRebin()
{
    Double_t* zero = new Double_t[cycles];
    Double_t* binNumArrDoubles = new Double_t [cycles];
    for(int i = 0; i < cycles; i++)
    {
        zero[i] = 0.0f;
    }
    for(int i = 0; i < cycles; i++)
    {
        binNumArrDoubles[i] = (Double_t) binNumArr[i];
    }

    for(int i = 0; i < numElements; i++)
    {
        fitDifferenceGraphs[i] = new TGraphErrors(cycles, binNumArrDoubles, meanDifferenceValues->GetHalfLifeArr(i), zero, meanDifferenceValues->GetHalfLifeErrorArr(i));
        fitDifferenceGraphs[i]->GetXaxis()->SetTitle("Number Bins");
        fitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit(s)");
        fitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Bateman Histo Mean Difference").c_str());

        singleFitDifferenceGraphs[i] = new TGraphErrors(cycles, binNumArrDoubles, singleMeanDifferenceValues->GetHalfLifeArr(i, i), zero, singleMeanDifferenceValues->GetHalfLifeErrorArr(i, i));
        singleFitDifferenceGraphs[i]->GetXaxis()->SetTitle("Number Bins");
        singleFitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit(s)");
        singleFitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Single Histo Mean Difference").c_str());
    }
    delete [] zero;
    delete [] binNumArrDoubles;
}

//creates and fills graphs for the difference in mean results for time difference
void Cycle::genMeanDifferenceGraphsTimeDifference()
{
    Double_t* zero = new Double_t[cycles];
    for(int i = 0; i < cycles; i++)
    {
        zero[i] = 0.0f;
    }

    for(int i = 0; i < numElements; i++)
    {
        fitDifferenceGraphs[i] = new TGraphErrors(cycles, timeArr, meanDifferenceValues->GetHalfLifeArr(i), zero, meanDifferenceValues->GetHalfLifeErrorArr(i));
        fitDifferenceGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        fitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit(s)");
        fitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Bateman Histo Mean Difference").c_str());

        singleFitDifferenceGraphs[i] = new TGraphErrors(cycles, timeArr, singleMeanDifferenceValues->GetHalfLifeArr(i, i), zero, singleMeanDifferenceValues->GetHalfLifeErrorArr(i, i));
        singleFitDifferenceGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        singleFitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit");
        singleFitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Single Histo Mean Difference(s)").c_str());
    }
    delete [] zero;
}

//calculates mean difference and error data for the single cycle data (dynamic)
void Cycle::genSingleMeanDifference()
{
    Double_t differenceHalfLife, differenceHalfLifeError;

    //calculates difference and error
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k < cycles; k++)
        {
            differenceHalfLife = batemanFitValues->GetAnHalfLife(k, i) - integralFitValues->GetAnHalfLife(k, i);
            meanDifferenceValues->SetAnHalfLife(k, i, differenceHalfLife);
            differenceHalfLife = singleBatemanFitValues->GetAnHalfLife(k, i, i) - singleIntegralFitValues->GetAnHalfLife(k, i, i);
            singleMeanDifferenceValues->SetAnHalfLife(k, i, i, differenceHalfLife);

            differenceHalfLifeError = sqrt(TMath::Power(batemanFitValues->GetAnHalfLifeError(k, i),2) + TMath::Power(integralFitValues->GetAnHalfLifeError(k, i),2));
            meanDifferenceValues->SetAnHalfLifeError(k, i, differenceHalfLifeError);
            differenceHalfLifeError = sqrt(TMath::Power(singleBatemanFitValues->GetAnHalfLifeError(k, i, i),2) + TMath::Power(singleIntegralFitValues->GetAnHalfLifeError(k, i, i),2));
            singleMeanDifferenceValues->SetAnHalfLifeError(k, i, i, differenceHalfLifeError);
        }
    }
}

//creates and fills the graphs for the seperate mean results for different bins
void Cycle::genSeperateMeanGraphsRebin()
{
    Double_t* tempFitVals, *tempFitErrors;
    Double_t* zero = new Double_t[cycles];
    Double_t* binNumArrDoubles = new Double_t [cycles];

    for(int i = 0; i < cycles; i++)
    {
        zero[i] = 0.0f;
    }
    for(int i = 0; i < cycles; i++)
    {
        binNumArrDoubles[i] = (Double_t) binNumArr[i];
    }

    for(int i = 0; i < numElements; i++)
    {
        tempFitVals = batemanFitValues->GetHalfLifeArr(i);
        tempFitErrors = batemanFitValues->GetHalfLifeErrorArr(i);
        batemanFitMeanGraphs[i] = new TGraphErrors(cycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        batemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        batemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        batemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Bateman Graph Mean").c_str());

        tempFitVals = singleBatemanFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleBatemanFitValues->GetHalfLifeErrorArr(i, i);
        singleBatemanFitMeanGraphs[i] = new TGraphErrors(cycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        singleBatemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        singleBatemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        singleBatemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single bateman Graph Mean").c_str());

        tempFitVals = integralFitValues->GetHalfLifeArr(i);
        tempFitErrors = integralFitValues->GetHalfLifeErrorArr(i);
        integralFitMeanGraphs[i] = new TGraphErrors(cycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        integralFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        integralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        integralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Integral Graph Mean").c_str());

        tempFitVals = singleIntegralFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleIntegralFitValues->GetHalfLifeErrorArr(i, i);
        singleIntegralFitMeanGraphs[i] = new TGraphErrors(cycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        singleIntegralFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        singleIntegralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(S)");
        singleIntegralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Integral Graph Mean").c_str());
    }

    delete [] zero;
    delete [] binNumArrDoubles;
}

//creates and fills the graphs for the seperate mean results (dynamic)
void Cycle::genSeperateMeanGraphsTimeDifference()
{
    Double_t* tempFitVals, *tempFitErrors;
    Double_t* zero = new Double_t[cycles];
    for(int i = 0; i < cycles; i++)
    {
        zero[i] = 0.0f;
    }

    for(int i = 0; i < numElements; i++)
    {
        tempFitVals = batemanFitValues->GetHalfLifeArr(i);
        tempFitErrors = batemanFitValues->GetHalfLifeErrorArr(i);
        batemanFitMeanGraphs[i] = new TGraphErrors(cycles, timeArr, tempFitVals, zero, tempFitErrors);
        batemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        batemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        batemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Bateman Graph Mean").c_str());

        tempFitVals = singleBatemanFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleBatemanFitValues->GetHalfLifeErrorArr(i, i);
        singleBatemanFitMeanGraphs[i] = new TGraphErrors(cycles, timeArr, tempFitVals, zero, tempFitErrors);
        singleBatemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        singleBatemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        singleBatemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Bateman Graph Mean").c_str());

        tempFitVals = integralFitValues->GetHalfLifeArr(i);
        tempFitErrors = integralFitValues->GetHalfLifeErrorArr(i);
        integralFitMeanGraphs[i] = new TGraphErrors(cycles, timeArr, tempFitVals, zero, tempFitErrors);
        integralFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        integralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        integralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Integral Graph Mean").c_str());

        tempFitVals = singleIntegralFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleIntegralFitValues->GetHalfLifeErrorArr(i, i);
        singleIntegralFitMeanGraphs[i] = new TGraphErrors(cycles, timeArr, tempFitVals, zero, tempFitErrors);
        singleIntegralFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        singleIntegralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        singleIntegralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Integral Graph Mean").c_str());
    }

    delete [] zero;
}

//runs the cycles, gets the mean for integral and bateman histograms and takes difference for the rebin function
void Cycle::runDifferenceMeanRebin()
{
    Double_t differenceHalfLife, differenceHalfLifeError;

    TH1D** multiRunHisto = decayChainRun->createRunResultHistos();
    TH1D** multiRunSingleHistos = decayChainRun->createRunResultHistosSingleElements();
    
    //elements in array stored in order, bateman, single
    for(int i = 0; i < cycles; i++)
    {
        //runs (numRuns) runs and generates run result histos
        decayChainRun->runNoChange(i);
        multiRunHisto = decayChainRun->fillRunResultHistos(multiRunHisto);
        multiRunSingleHistos = decayChainRun->fillRunResultHistosSingleElement(multiRunSingleHistos);

        //does mean difference calculation for each element in decay chain
        for(int k = 0; k < numElements; k++)
        {
            differenceHalfLife = multiRunHisto[(k*2)]->GetMean() - multiRunHisto[(k*2)+1]->GetMean();
            meanDifferenceValues->SetAnHalfLife(i, k, differenceHalfLife);
            differenceHalfLife = multiRunSingleHistos[(k*2)]->GetMean() - multiRunSingleHistos[(k*2)+1]->GetMean();
            singleMeanDifferenceValues->SetAnHalfLife(i, k, k, differenceHalfLife);

            differenceHalfLifeError = sqrt(pow(multiRunHisto[(k*2)+1]->GetMeanError(),2) + pow(multiRunHisto[(k*2)]->GetMeanError(),2));
            meanDifferenceValues->SetAnHalfLifeError(i, k, differenceHalfLifeError);
            differenceHalfLifeError = sqrt(pow(multiRunSingleHistos[(k*2)+1]->GetMeanError(),2) + pow(multiRunSingleHistos[(k*2)]->GetMeanError(),2));
            singleMeanDifferenceValues->SetAnHalfLifeError(i, k, k, differenceHalfLifeError);
        }
    }

    for(int i = 0; i < numElements*2; i++)
    {
        delete multiRunHisto[i];
        delete multiRunSingleHistos[i];
    }
    delete [] multiRunHisto;
    delete [] multiRunSingleHistos;
}

//runs the cycles, gets the mean for integral and bateman histograms and takes difference (dynamic)
void Cycle::runDifferenceMeanTimeDifference()
{
    Double_t differenceHalfLife, differenceHalfLifeError;

    TH1D** multiRunHisto = decayChainRun->createRunResultHistos();
    TH1D** multiRunSingleHistos = decayChainRun->createRunResultHistosSingleElements();
    
    //elements in array stored in order, bateman, single
    for(int i = 0; i < cycles; i++)
    {
        //runs (numRuns) runs and generates run result histos
        decayChainRun->runNoChange(i);
        multiRunHisto = decayChainRun->fillRunResultHistos(multiRunHisto);
        multiRunSingleHistos = decayChainRun->fillRunResultHistosSingleElement(multiRunSingleHistos);

        //does mean difference calculation for each element in decay chain
        for(int k = 0; k < numElements; k++)
        {
            differenceHalfLife = multiRunHisto[(k*2)]->GetMean() - multiRunHisto[(k*2)+1]->GetMean();
            meanDifferenceValues->SetAnHalfLife(i, k, differenceHalfLife);
            differenceHalfLife = multiRunSingleHistos[(k*2)]->GetMean() - multiRunSingleHistos[(k*2)+1]->GetMean();
            singleMeanDifferenceValues->SetAnHalfLife(i, k, k, differenceHalfLife);

            differenceHalfLifeError = sqrt(pow(multiRunHisto[(k*2)+1]->GetMeanError(),2) + pow(multiRunHisto[(k*2)]->GetMeanError(),2));
            meanDifferenceValues->SetAnHalfLifeError(i, k, differenceHalfLifeError);
            differenceHalfLifeError = sqrt(pow(multiRunSingleHistos[(k*2)+1]->GetMeanError(),2) + pow(multiRunSingleHistos[(k*2)]->GetMeanError(),2));
            singleMeanDifferenceValues->SetAnHalfLifeError(i, k, k, differenceHalfLifeError);
        }
    }

    for(int i = 0; i < numElements*2; i++)
    {
        delete multiRunHisto[i];
        delete multiRunSingleHistos[i];
    }
    delete [] multiRunHisto;
    delete [] multiRunSingleHistos;
}

//runs the cycles and puts the data of the integral method mean and the bateman method mean into their respective arrays (dynamic)
void Cycle::runSeperateMeanRebin()
{
    Double_t tempMeanVal;
    Double_t tempErrorVal;

    TH1D** multiRunHisto = decayChainRun->createRunResultHistos();
    TH1D** multiRunHistoSingle = decayChainRun->createRunResultHistosSingleElements();

    //running the cycle and putting histograms means in respective arrays
    for(int i = 0; i < cycles; i++)
    {
        decayChainRun->runNoChange(i);
        //fills histograms with results from multiple runs in order to extract means
        multiRunHisto = decayChainRun->fillRunResultHistos(multiRunHisto);
        multiRunHistoSingle = decayChainRun->fillRunResultHistosSingleElement(multiRunHistoSingle);

        //moves histo means into respective data storages
        for(int k = 0; k < numElements; k++)
        {
            tempMeanVal = multiRunHisto[(k*2)]->GetMean();
            batemanFitValues->SetAnHalfLife(i, k, tempMeanVal);
            tempErrorVal = multiRunHisto[(k*2)]->GetMeanError();
            batemanFitValues->SetAnHalfLifeError(i, k, tempErrorVal);
            tempMeanVal = multiRunHisto[(k*2)+1]->GetMean();
            integralFitValues->SetAnHalfLife(i, k, tempMeanVal);
            tempErrorVal = multiRunHisto[(k*2)+1]->GetMeanError();
            integralFitValues->SetAnHalfLifeError(i, k, tempErrorVal);

            tempMeanVal = multiRunHistoSingle[(k*2)]->GetMean();
            singleBatemanFitValues->SetAnHalfLife(i, k, k, tempMeanVal);
            tempErrorVal = multiRunHistoSingle[(k*2)]->GetMeanError();
            singleBatemanFitValues->SetAnHalfLifeError(i, k, k, tempErrorVal);
            tempMeanVal = multiRunHistoSingle[(k*2)+1]->GetMean();
            singleIntegralFitValues->SetAnHalfLife(i, k, k, tempMeanVal);
            tempErrorVal = multiRunHistoSingle[(k*2)+1]->GetMeanError();
            singleIntegralFitValues->SetAnHalfLifeError(i, k, k, tempErrorVal);
        }
    }

    for(int i = 0; i < numElements*2; i++)
    {
        delete multiRunHisto[i];
        delete multiRunHistoSingle[i];
    }
    delete [] multiRunHisto;
    delete [] multiRunHistoSingle;
}

//runs the cycles and puts the data of the integral method mean and the bateman method mean into their respective arrays (dynamic)
void Cycle::runSeperateMeanTimeDifference()
{
    Double_t tempMeanVal, tempErrorVal;
    Double_t *runningSingleBateman = new Double_t [numElements];
    Double_t *runningTotalBateman = new Double_t [numElements];
    Double_t *runningSingleIntegral = new Double_t [numElements];
    Double_t *runningTotalIntegral = new Double_t [numElements];
    for(int i = 0; i < numElements; i++)
    {
        runningTotalIntegral[i] = 0.0f;
        runningTotalBateman[i] = 0.0f;
        runningSingleBateman[i] = 0.0f;
        runningSingleIntegral[i] = 0.0f;
    }

    TH1D** multiRunHisto = decayChainRun->createRunResultHistos();
    TH1D** multiRunHistoSingle = decayChainRun->createRunResultHistosSingleElements();

    //running the cycle and putting histograms means in respective arrays
    for(int i = 0; i < cycles; i++)
    {
        decayChainRun->runNoChange(i);
        //fills histograms with results from multiple runs in order to extract means
        multiRunHisto = decayChainRun->fillRunResultHistos(multiRunHisto);
        multiRunHistoSingle = decayChainRun->fillRunResultHistosSingleElement(multiRunHistoSingle);

        //moves histo means into respective data storages
        for(int k = 0; k < numElements; k++)
        {
            tempMeanVal = multiRunHisto[(k*2)]->GetMean();
            batemanFitValues->SetAnHalfLife(i, k, tempMeanVal);
            runningTotalBateman[k] += tempMeanVal;
            tempErrorVal = multiRunHisto[(k*2)]->GetMeanError();
            batemanFitValues->SetAnHalfLifeError(i, k, tempErrorVal);
            tempMeanVal = multiRunHisto[(k*2)+1]->GetMean();
            integralFitValues->SetAnHalfLife(i, k, tempMeanVal);
            runningTotalIntegral[k] += tempMeanVal;
            tempErrorVal = multiRunHisto[(k*2)+1]->GetMeanError();
            integralFitValues->SetAnHalfLifeError(i, k, tempErrorVal);

            tempMeanVal = multiRunHistoSingle[(k*2)]->GetMean();
            singleBatemanFitValues->SetAnHalfLife(i, k, k, tempMeanVal);
            runningSingleBateman[k] += tempMeanVal;
            tempErrorVal = multiRunHistoSingle[(k*2)]->GetMeanError();
            singleBatemanFitValues->SetAnHalfLifeError(i, k, k, tempErrorVal);
            tempMeanVal = multiRunHistoSingle[(k*2)+1]->GetMean();
            runningSingleIntegral[k] += tempMeanVal;
            singleIntegralFitValues->SetAnHalfLife(i, k, k, tempMeanVal);
            tempErrorVal = multiRunHistoSingle[(k*2)+1]->GetMeanError();
            singleIntegralFitValues->SetAnHalfLifeError(i, k, k, tempErrorVal);
        }
    }

    for(int i = 0; i < numElements; i++)
    {
        cout << elementStrNames[i] << " AVERAGE BATEMAN TOTAL: " << runningTotalBateman[i] / cycles << endl;
        cout << elementStrNames[i] << " AVERAGE INTEGRAL TOTAL: " << runningTotalIntegral[i] / cycles << endl;
        cout << elementStrNames[i] << " AVERAGE BATEMAN SINGLE: " << runningSingleBateman[i] / cycles << endl;
        cout << elementStrNames[i] << " AVERAGE INTEGRAL SINGLE: " << runningSingleIntegral[i] / cycles << endl << endl;
    }

    delete [] runningTotalBateman;
    delete [] runningTotalIntegral;
    delete [] runningSingleBateman;
    delete [] runningSingleIntegral;
    for(int i = 0; i < numElements*2; i++)
    {
        delete multiRunHisto[i];
        delete multiRunHistoSingle[i];
    }
    delete [] multiRunHisto;
    delete [] multiRunHistoSingle;
}


//does runs cycle for the single histogram (dynamic)
void Cycle::runSeperateSingleGen()
{
    Double_t N0;
    Double_t N0Error;
    Double_t halfLife;
    Double_t halfLifeError;

    ChainRunFitValues* tempVals;
    SingleChainRunFitValues* singleTempVals;

    //running the cycle and puts the fit values in the arrays
    for(int i = 0; i < cycles; i++)
    {
        //runs the fit
        decayChainRun->runNoChangeGenOnce(i, 0);
        
        //gets total bateman fit data from fit
        tempVals = decayChainRun->getBatemanFitValues();

        //move total bateman parameters in to respective class
        for(int k = 0; k < numElements; k++)
        {
            N0 = tempVals->GetAnN0(0, k);
            N0Error = tempVals->GetAnN0Error(0, k);
            halfLife = tempVals->GetAnHalfLife(0, k);
            halfLifeError = tempVals->GetAnHalfLifeError(0, k);

            batemanFitValues->SetAnN0(i, k, N0);
            batemanFitValues->SetAnN0Error(i, k, N0Error);
            batemanFitValues->SetAnHalfLife(i, k, halfLife);
            batemanFitValues->SetAnHalfLifeError(i, k, halfLifeError);
        }
        //gets integral fit data from fit
        tempVals = decayChainRun->getIntegralFitValues();

        //move total integral parameters in to respective class
        for(int k = 0; k < numElements; k++)
        {
            N0 = tempVals->GetAnN0(0, k);
            N0Error = tempVals->GetAnN0Error(0, k);
            halfLife = tempVals->GetAnHalfLife(0, k);
            halfLifeError = tempVals->GetAnHalfLifeError(0, k);

            integralFitValues->SetAnN0Error(i, k, N0Error);
            integralFitValues->SetAnHalfLife(i, k, halfLife);
            integralFitValues->SetAnN0(i, k, N0);
            integralFitValues->SetAnHalfLifeError(i, k, halfLifeError);
        }

        //gets single bateman fit data from fit
        singleTempVals = decayChainRun->getSingleBatemanFitValues();

        //move single bateman parameters in to respective class
        for(int k = 0; k < numElements; k++)
        {
            for(int subElement = 0; subElement < k+1; subElement++)
            {
                N0 = singleTempVals->GetAnN0(0, k, subElement);
                N0Error = singleTempVals->GetAnN0Error(0, k, subElement);
                halfLife = singleTempVals->GetAnHalfLife(0, k, subElement);
                halfLifeError = singleTempVals->GetAnHalfLifeError(0, k, subElement);

                singleBatemanFitValues->SetAnN0(i, k, subElement, N0);
                singleBatemanFitValues->SetAnN0Error(i, k, subElement, N0Error);
                singleBatemanFitValues->SetAnHalfLife(i, k, subElement, halfLife);
                singleBatemanFitValues->SetAnHalfLifeError(i, k, subElement, halfLifeError);
            }
        }

        //gets single bateman fit data from fit
        singleTempVals = decayChainRun->getSingleIntegralFitValues();

        //move single bateman parameters in to respective class
        for(int k = 0; k < numElements; k++)
        {
            for(int subElement = 0; subElement < k+1; subElement++)
            {
                N0 = singleTempVals->GetAnN0(0, k, subElement);
                N0Error = singleTempVals->GetAnN0Error(0, k, subElement);
                halfLife = singleTempVals->GetAnHalfLife(0, k, subElement);
                halfLifeError = singleTempVals->GetAnHalfLifeError(0, k, subElement);

                singleIntegralFitValues->SetAnN0(i, k, subElement, N0);
                singleIntegralFitValues->SetAnN0Error(i, k, subElement, N0Error);
                singleIntegralFitValues->SetAnHalfLife(i, k, subElement, halfLife);
                singleIntegralFitValues->SetAnHalfLifeError(i, k, subElement, halfLifeError);
            }
        }
    }
}

#endif