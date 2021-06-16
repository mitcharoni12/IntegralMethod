#ifndef CYCLE_H
#define CYCLE_H

#include <iostream>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Run.h"
#include "ElementFit.h"
#include "RunFitValues.h"
#include "ChainRunFitValues.h"
#include "SingleChainRunFitValues.h"

using namespace std;

class Cycle{
    private:
        Int_t cycles, numElements, incChoice;
        Double_t x_inc, x_start, x_stop;
        Double_t* timeArr;
        Run* decayChainRun;
        ElementFit* element;
        string* elementStrNames;
        ChainRunFitValues* regularFitValues, *integralFitValues;
        SingleChainRunFitValues* singleRegularFitValues, *singleIntegralFitValues;
        TGraphErrors** regularFitMeanGraphs, **integralFitMeanGraphs, **singleRegularFitMeanGraphs
        ,**singleIntegralFitMeanGraphs, **fitDifferenceGraphs, **singleFitDifferenceGraphs;
        //TCanvas* test = new TCanvas("test", "test", 500, 500);
        //TCanvas* test2 = new TCanvas("test2", "test2", 500, 500);
        //helper functions
    public:
        Cycle(Int_t cycles, Run* decayChainRun, ElementFit* element, Double_t x_start, Double_t x_stop, Double_t x_inc, Int_t incChoice);
        ~Cycle();
        void displayMeanDifferenceGraphs(TCanvas** canvasArr, TGraphErrors** cycleMeanDifferenceGraphs);
        void displayMeanSeperateGraphs(TCanvas** canvasArr);
        //TGraphErrors** genMeanDifferenceGraphs(resultStorage<Double_t>** cycleMeanDifference);
        void genSeperateMeanGraphs();
        //resultStorage<Double_t>** genSingleMeanDifference(resultStorage<Double_t>** cycleSingleGenResult);
        //resultStorage<Double_t>** runDifferenceMean();
        void runSeperateMean();
        void runSeperateSingleGen();
};

Cycle::Cycle(Int_t cycles, Run* decayChainRun, ElementFit* element, Double_t x_start, Double_t x_stop, Double_t x_inc, Int_t incChoice)
{
    this->cycles = cycles;
    this->element = element;
    this->decayChainRun = decayChainRun;
    timeArr = new Double_t [cycles];
    numElements = element->getNumElements();
    this->x_start = x_start;
    this->x_stop = x_stop;
    this->x_inc = x_inc;
    this->incChoice = incChoice;
    element->setNumRuns(decayChainRun->getNumRuns());
    element->setNumCycles(cycles);
    elementStrNames = decayChainRun->getElementStringNames();
    regularFitValues = new ChainRunFitValues(numElements, cycles);
    integralFitValues = new ChainRunFitValues(numElements, cycles);
    singleRegularFitValues = new SingleChainRunFitValues(numElements, cycles);
    singleIntegralFitValues = new SingleChainRunFitValues(numElements, cycles);
    regularFitMeanGraphs = new TGraphErrors* [numElements];
    integralFitMeanGraphs = new TGraphErrors* [numElements];
    singleRegularFitMeanGraphs = new TGraphErrors* [numElements];
    singleIntegralFitMeanGraphs = new TGraphErrors* [numElements];
    fitDifferenceGraphs = new TGraphErrors* [numElements];
    singleFitDifferenceGraphs = new TGraphErrors* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        regularFitMeanGraphs[i] = nullptr;
        integralFitMeanGraphs[i] = nullptr;
        singleRegularFitMeanGraphs[i] = nullptr;
        singleIntegralFitMeanGraphs[i] = nullptr;
        fitDifferenceGraphs[i] = nullptr;
        singleFitDifferenceGraphs[i] = nullptr;
    }
}

Cycle::~Cycle()
{
    delete [] timeArr;
    delete [] regularFitMeanGraphs;
    delete [] integralFitMeanGraphs;
    delete [] singleRegularFitMeanGraphs;
    delete [] singleIntegralFitMeanGraphs;
    delete [] fitDifferenceGraphs;
    delete [] singleFitDifferenceGraphs;
    delete regularFitValues;
    delete integralFitValues;
    delete singleRegularFitValues;
    delete singleIntegralFitValues;
}

/*
//displays the graphs for the mean difference
void Cycle::displayMeanDifferenceGraphs(TCanvas** canvasArr, TGraphErrors** cycleMeanDifferenceGraphs)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(1);
        cycleMeanDifferenceGraphs[i*2]->Draw();
        canvasArr[i]->cd(2);
        cycleMeanDifferenceGraphs[(i*2)+1]->Draw();
    }
}
*/

//displays the graphs for the seperate mean
void Cycle::displayMeanSeperateGraphs(TCanvas** canvasArr)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(1);
        regularFitMeanGraphs[i]->Draw();
        canvasArr[i]->cd(2);
        integralFitMeanGraphs[i]->Draw();
        canvasArr[i]->cd(3);
        singleRegularFitMeanGraphs[i]->Draw();
        canvasArr[i]->cd(4);
        singleIntegralFitMeanGraphs[i]->Draw();
    }
}

/*
//creates and fills graphs for the difference in mean results (dynamic)
TGraphErrors** Cycle::genMeanDifferenceGraphs(resultStorage<Double_t>** cycleMeanDifference)
{
    Double_t* zero = new Double_t[cycles];
    for(int i = 0; i < cycles; i++)
    {
        zero[i] = 0.0f;
    }

    int numRegPlusSingleHistos = numElements*2;
    TGraphErrors** cycleMeanDifferenceGraphs = new TGraphErrors* [numRegPlusSingleHistos];
    for(int i = 0; i < numElements; i++)
    {
        cycleMeanDifferenceGraphs[(i*2)] = new TGraphErrors(cycles, timeArr, (cycleMeanDifference[0]->getDoubleArrStorage())[i], zero, (cycleMeanDifference[2]->getDoubleArrStorage())[i]);
        cycleMeanDifferenceGraphs[(i*2)]->GetXaxis()->SetTitle("Time");
        cycleMeanDifferenceGraphs[(i*2)]->GetYaxis()->SetTitle("IntegralFit - RegularFit");
        cycleMeanDifferenceGraphs[(i*2)]->SetTitle((elementStrNames[i] + " Regular Histo Mean Difference").c_str());

        cycleMeanDifferenceGraphs[(i*2)+1] = new TGraphErrors(cycles, timeArr, (cycleMeanDifference[1]->getDoubleArrStorage())[i], zero, (cycleMeanDifference[3]->getDoubleArrStorage())[i]);
        cycleMeanDifferenceGraphs[(i*2)+1]->GetXaxis()->SetTitle("Time");
        cycleMeanDifferenceGraphs[(i*2)+1]->GetYaxis()->SetTitle("IntegralFit - RegularFit");
        cycleMeanDifferenceGraphs[(i*2)+1]->SetTitle((elementStrNames[i] + " Single Histo Mean Difference").c_str());
    }
    delete [] zero;
    return cycleMeanDifferenceGraphs;
}
*/

/*
//calculates mean difference and error data for the single cycle data (dynamic)
resultStorage<Double_t>** Cycle::genSingleMeanDifference(resultStorage<Double_t>** cycleSingleGenResult)
{
    resultStorage<Double_t>** meanDifferenceStore = new resultStorage<Double_t>* [4];
    meanDifferenceStore[0] = new resultStorage<Double_t>(2, cycles, numElements); //total mean difference
    meanDifferenceStore[1] = new resultStorage<Double_t>(2, cycles, numElements); //single mean difference
    meanDifferenceStore[2] = new resultStorage<Double_t>(2, cycles, numElements); //total mean differenc error
    meanDifferenceStore[3] = new resultStorage<Double_t>(2, cycles, numElements); //single mean difference error
    //calculates difference and error
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k < cycles; k++)
        {
            (meanDifferenceStore[0]->getDoubleArrStorage())[i][k] = (cycleSingleGenResult[2]->getDoubleArrStorage())[i][k] - (cycleSingleGenResult[0]->getDoubleArrStorage())[i][k]; 
            (meanDifferenceStore[2]->getDoubleArrStorage())[i][k] = (cycleSingleGenResult[6]->getDoubleArrStorage())[i][k] - (cycleSingleGenResult[4]->getDoubleArrStorage())[i][k];
            (meanDifferenceStore[1]->getDoubleArrStorage())[i][k] = sqrt(pow((cycleSingleGenResult[3]->getDoubleArrStorage())[i][k],2) + pow((cycleSingleGenResult[1]->getDoubleArrStorage())[i][k],2));
            (meanDifferenceStore[3]->getDoubleArrStorage())[i][k] = sqrt(pow((cycleSingleGenResult[7]->getDoubleArrStorage())[i][k],2) + pow((cycleSingleGenResult[5]->getDoubleArrStorage())[i][k],2));
        }
    }
    return meanDifferenceStore;
}
*/

//creates and fills the graphs for the seperate mean results (dynamic)
void Cycle::genSeperateMeanGraphs()
{
    Double_t* tempFitVals, *tempFitErrors;
    Double_t* zero = new Double_t[cycles];
    for(int i = 0; i < cycles; i++)
    {
        zero[i] = 0.0f;
    }
    TGraphErrors** meanSeperateGraphs = new TGraphErrors* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        tempFitVals = regularFitValues->GetHalfLifeArr(i);
        tempFitErrors = regularFitValues->GetHalfLifeErrorArr(i);
        regularFitMeanGraphs[i] = new TGraphErrors(cycles, timeArr, tempFitVals, zero, tempFitErrors);
        regularFitMeanGraphs[i]->GetXaxis()->SetTitle("Time");
        regularFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value");
        regularFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Regular Graph Mean").c_str());

        tempFitVals = singleRegularFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleRegularFitValues->GetHalfLifeErrorArr(i, i);
        singleRegularFitMeanGraphs[i] = new TGraphErrors(cycles, timeArr, tempFitVals, zero, tempFitErrors);
        singleRegularFitMeanGraphs[i]->GetXaxis()->SetTitle("Time");
        singleRegularFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value");
        singleRegularFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Regular Graph Mean").c_str());

        tempFitVals = integralFitValues->GetHalfLifeArr(i);
        tempFitErrors = integralFitValues->GetHalfLifeErrorArr(i);
        integralFitMeanGraphs[i] = new TGraphErrors(cycles, timeArr, tempFitVals, zero, tempFitErrors);
        integralFitMeanGraphs[i]->GetXaxis()->SetTitle("Time");
        integralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value");
        integralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Integral Graph Mean").c_str());

        tempFitVals = singleIntegralFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleIntegralFitValues->GetHalfLifeErrorArr(i, i);
        singleIntegralFitMeanGraphs[i] = new TGraphErrors(cycles, timeArr, tempFitVals, zero, tempFitErrors);
        singleIntegralFitMeanGraphs[i]->GetXaxis()->SetTitle("Time");
        singleIntegralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value");
        singleIntegralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Integral Graph Mean").c_str());

    }

    delete [] zero;
}

    /*
//runs the cycles, gets the mean for integral and regular histograms and takes difference (dynamic)
resultStorage<Double_t>** Cycle::runDifferenceMean()
{
    //setting inital values for time, 1= start const end moving, 2= end const start moving
    element->setTimeRunStart(x_start);
    element->setTimeRunEnd(x_stop);

    resultStorage<Double_t>** meanDifferenceData = new resultStorage<Double_t>* [4];
    meanDifferenceData[0] = new resultStorage<Double_t>(2, cycles, numElements);//Total Mean Difference
    meanDifferenceData[1] = new resultStorage<Double_t>(2, cycles, numElements);//Single Mean Difference
    meanDifferenceData[2] = new resultStorage<Double_t>(2, cycles, numElements);//Total Mean Difference Error
    meanDifferenceData[3] = new resultStorage<Double_t>(2, cycles, numElements);//Single Mean Difference Error

    TH1D** multiRunHisto = decayChainRun->createRunResultHistos();
    TH1D** multiRunSingleHistos = decayChainRun->createRunResultHistosSingleElements();
    
    //elements in array stored in order, regular, single
    // 0 = total mean difference 1 = single mean difference 2 = total mean difference error 3 = single mean difference error
    for(int i = 0; i < cycles; i++)
    {
        //runs (numRuns) runs and generates run result histos
        decayChainRun->runNoChange();
        multiRunHisto = decayChainRun->fillRunResultHistos(multiRunHisto);
        multiRunSingleHistos = decayChainRun->fillRunResultHistosSingleElement(multiRunSingleHistos);
        //does mean difference calculation for each element in decay chain
        for(int k = 0; k < numElements; k++)
        {
            (meanDifferenceData[0]->getDoubleArrStorage())[k][i] = multiRunHisto[(k*2)+1]->GetMean() - multiRunHisto[(k*2)]->GetMean();
            (meanDifferenceData[1]->getDoubleArrStorage())[k][i] = multiRunSingleHistos[(k*2)+1]->GetMean() - multiRunSingleHistos[(k*2)]->GetMean();
            (meanDifferenceData[2]->getDoubleArrStorage())[k][i] = sqrt(pow(multiRunHisto[(k*2)+1]->GetMeanError(),2) - pow(multiRunHisto[(k*2)]->GetMeanError(),2));
            (meanDifferenceData[3]->getDoubleArrStorage())[k][i] = sqrt(pow(multiRunSingleHistos[(k*2)+1]->GetMeanError(),2) - pow(multiRunSingleHistos[(k*2)]->GetMeanError(),2));
        }
        //sets next range for time and sets current time in array, 1= start const end moving, 2= end const start moving
        if(incChoice == 1)
        {
            element->setTimeRunEnd((i+1.0)*x_inc + x_stop);
            timeArr[i] = ((i*x_inc) + x_stop);
        }else if(incChoice == 2)
        {
            element->setTimeRunStart((i+1.0)*x_inc + x_start);
            timeArr[i] = ((i*x_inc) + x_start);
        }
    }

    return meanDifferenceData;
    return nullptr;
}
    */

//runs the cycles and puts the data of the integral method mean and the regular method mean into their respective arrays (dynamic)
void Cycle::runSeperateMean()
{
    //setting inital values for time, 1= start const end moving, 2= end const start moving
    element->setTimeRunStart(x_start);
    element->setTimeRunEnd(x_stop);
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
            regularFitValues->SetAnHalfLifeError(i, k, tempMeanVal);
            tempErrorVal = multiRunHisto[(k*2)]->GetMeanError();
            regularFitValues->SetAnHalfLife(i, k, tempErrorVal);
            tempMeanVal = multiRunHisto[(k*2)+1]->GetMean();
            integralFitValues->SetAnHalfLifeError(i, k, tempMeanVal);
            tempErrorVal = multiRunHisto[(k*2)+1]->GetMeanError();
            integralFitValues->SetAnHalfLife(i, k, tempErrorVal);

            tempMeanVal = multiRunHistoSingle[(k*2)]->GetMean();
            singleRegularFitValues->SetAnHalfLife(i, k, k, tempMeanVal);
            tempErrorVal = multiRunHistoSingle[(k*2)]->GetMean();
            singleRegularFitValues->SetAnHalfLifeError(i, k, k, tempErrorVal);
            tempMeanVal = multiRunHistoSingle[(k*2)+1]->GetMean();
            singleIntegralFitValues->SetAnHalfLife(i, k, k, tempMeanVal);
            tempErrorVal = multiRunHistoSingle[(k*2)+1]->GetMean();
            singleIntegralFitValues->SetAnHalfLifeError(i, k, k, tempErrorVal);
        }

        //sets next range for time and sets current time in array, 1= start const end moving, 2= end const start moving
        if(incChoice == 1)
        {
            element->setTimeRunEnd((i+1.0)*x_inc + x_stop);
            timeArr[i] = ((i*x_inc) + x_stop);
        }else if(incChoice == 2)
        {
            element->setTimeRunStart((i+1.0)*x_inc + x_start);
            timeArr[i] = ((i*x_inc) + x_start);
        }
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

    //setting inital values for time, 1 = start const end moving, 2 = start moving end const
    element->setTimeRunStart(x_start);
    element->setTimeRunEnd(x_stop);

    ChainRunFitValues* tempVals;
    SingleChainRunFitValues* singleTempVals;

    //running the cycle and puts the fit values in the arrays
    for(int i = 0; i < cycles; i++)
    {
        //runs the fit
        decayChainRun->runNoChangeGenOnce(i, 0);
        
        //gets total regular fit data from fit
        tempVals = decayChainRun->getRegularFitValues();

        //move total regular parameters in to respective class
        for(int k = 0; k < numElements; k++)
        {
            N0 = tempVals->GetAnN0(0, k);
            N0Error = tempVals->GetAnN0Error(0, k);
            halfLife = tempVals->GetAnHalfLife(0, k);
            halfLifeError = tempVals->GetAnHalfLifeError(0, k);

            regularFitValues->SetAnN0(i, k, N0);
            regularFitValues->SetAnN0Error(i, k, N0Error);
            regularFitValues->SetAnHalfLife(i, k, halfLife);
            regularFitValues->SetAnHalfLifeError(i, k, halfLifeError);
        }
        //gets integral fit data from fit
        tempVals = decayChainRun->getRegularFitValues();

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

        //gets single regular fit data from fit
        singleTempVals = decayChainRun->getSingleRegularFitValues();

        //move single regular parameters in to respective class
        for(int k = 0; k < numElements; k++)
        {
            for(int subElement = 0; subElement < k+1; subElement++)
            {
                N0 = singleTempVals->GetAnN0(0, k, subElement);
                N0Error = singleTempVals->GetAnN0Error(0, k, subElement);
                halfLife = singleTempVals->GetAnHalfLife(0, k, subElement);
                halfLifeError = singleTempVals->GetAnHalfLifeError(0, k, subElement);

                singleRegularFitValues->SetAnN0(i, k, subElement, N0);
                singleRegularFitValues->SetAnN0Error(i, k, subElement, N0Error);
                singleRegularFitValues->SetAnHalfLife(i, k, subElement, halfLife);
                singleRegularFitValues->SetAnHalfLifeError(i, k, subElement, halfLifeError);
            }
        }

        //gets single regular fit data from fit
        singleTempVals = decayChainRun->getSingleIntegralFitValues();

        //move single regular parameters in to respective class
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

        //sets next range for time and sets current time in array, 1= start const end moving, 2= end const start moving
        if(incChoice == 1)
        {
            element->setTimeRunEnd((i+1.0)*x_inc + x_stop);
            timeArr[i] = ((i*x_inc) + x_stop);
        }else if(incChoice == 2)
        {
            element->setTimeRunStart((i+1.0)*x_inc + x_start);
            timeArr[i] = ((i*x_inc) + x_start);
        }
    }
}

#endif