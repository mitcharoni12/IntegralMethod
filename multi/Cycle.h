#ifndef CYCLE_H
#define CYCLE_H

#include <iostream>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Run.h"
#include "ElementFit.h"
#include "resultStorage.h"
#include "RunFitValues.h"
#include "ChainRunFitValues.h"
#include "SingleChainRunFitVales.h"

using namespace std;

class Cycle{
    private:
        Int_t cycles, numElements, incChoice;
        Double_t x_inc, x_start, x_stop;
        Double_t* timeArr;
        Run* decayChainRun;
        ElementFit* element;
        string* elementStrNames;
        ChainRunFitValues* regularFitValues;
        ChainRunFitValues* integralFitValues;
        SingleChainRunFitVales* singleRegularFitValues;
        SingleChainRunFitVales* singleIntegralFitValues;
        TCanvas* testCan;
        //helper functions
    public:
        Cycle(Int_t cycles, Run* decayChainRun, ElementFit* element, Double_t x_start, Double_t x_stop, Double_t x_inc, Int_t incChoice);
        ~Cycle();
        void displayMeanDifferenceGraphs(TCanvas** canvasArr, TGraphErrors** cycleMeanDifferenceGraphs);
        void displayMeanSeperateGraphs(TCanvas** canvasArr, TGraphErrors** meanResultGraphs);
        TGraphErrors** genMeanDifferenceGraphs(resultStorage<Double_t>** cycleMeanDifference);
        TGraphErrors** genSeperateMeanGraphs(resultStorage<Double_t>** cycleMeanResults);
        resultStorage<Double_t>** genSingleMeanDifference(resultStorage<Double_t>** cycleSingleGenResult);
        resultStorage<Double_t>** runDifferenceMean();
        resultStorage<Double_t>** runSeperateMean();
        resultStorage<Double_t>** runSeperateSingleGen();
};

Cycle::Cycle(Int_t cycles, Run* decayChainRun, ElementFit* element, Double_t x_start, Double_t x_stop, Double_t x_inc, Int_t incChoice)
{
    this->cycles = cycles;
    this->element = element;
    this->decayChainRun = decayChainRun;
    timeArr = new Double_t [cycles];
    numElements = element->getNumElements();
    testCan = new TCanvas("testCan", "testCan", 500, 500);
    this->x_start = x_start;
    this->x_stop = x_stop;
    this->x_inc = x_inc;
    this->incChoice = incChoice;
    elementStrNames = decayChainRun->getElementStringNames();
    /*
    regularFitValues = new ChainRunFitValues(numElements, cycles);
    integralFitValues = new ChainRunFitValues(numElements, cycles);
    singleRegularFitValues = new SingleChainRunFitValues(numElements, cycles);
    singleIntegralFitValues = new SingleChainRunFitValues(numElements, cycles);
    */

}

Cycle::~Cycle()
{
    delete [] timeArr;
}

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

//displays the graphs for the seperate mean
void Cycle::displayMeanSeperateGraphs(TCanvas** canvasArr, TGraphErrors** meanResultGraphs)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(1);
        meanResultGraphs[(i*4)]->Draw();
        canvasArr[i]->cd(2);
        meanResultGraphs[(i*4)+1]->Draw();
        canvasArr[i]->cd(3);
        meanResultGraphs[(i*4)+2]->Draw();
        canvasArr[i]->cd(4);
        meanResultGraphs[(i*4)+3]->Draw();
    }
}

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

//creates and fills the graphs for the seperate mean results (dynamic)
TGraphErrors** Cycle::genSeperateMeanGraphs(resultStorage<Double_t>** cycleMeanResults)
{
    Double_t* zero = new Double_t[cycles];
    for(int i = 0; i < cycles; i++)
    {
        zero[i] = 0.0f;
    }
    TGraphErrors** meanSeperateGraphs = new TGraphErrors* [numElements];
    //0 = regular graph 1 = single regular graph 2 = integral graph 3 = single integral graph
    for(int i = 0; i < numElements; i++)
    {
        meanSeperateGraphs[(i*4)] = new TGraphErrors(cycles, timeArr, (cycleMeanResults[0]->getDoubleArrStorage())[i], zero, (cycleMeanResults[1]->getDoubleArrStorage())[i]);
        meanSeperateGraphs[(i*4)]->GetXaxis()->SetTitle("Time");
        meanSeperateGraphs[(i*4)]->GetYaxis()->SetTitle("Fit Value");
        meanSeperateGraphs[(i*4)]->SetTitle((elementStrNames[i] + " Regular Graph Mean").c_str());

        meanSeperateGraphs[(i*4)+1] = new TGraphErrors(cycles, timeArr, (cycleMeanResults[4]->getDoubleArrStorage())[i], zero, (cycleMeanResults[5]->getDoubleArrStorage())[i]);
        meanSeperateGraphs[(i*4)+1]->GetXaxis()->SetTitle("Time");
        meanSeperateGraphs[(i*4)+1]->GetYaxis()->SetTitle("Fit Value");
        meanSeperateGraphs[(i*4)+1]->SetTitle((elementStrNames[i] + " Single Regular Graph Mean").c_str());

        meanSeperateGraphs[(i*4)+2] = new TGraphErrors(cycles, timeArr, (cycleMeanResults[2]->getDoubleArrStorage())[i], zero, (cycleMeanResults[3]->getDoubleArrStorage())[i]);
        meanSeperateGraphs[(i*4)+2]->GetXaxis()->SetTitle("Time");
        meanSeperateGraphs[(i*4)+2]->GetYaxis()->SetTitle("Fit Value");
        meanSeperateGraphs[(i*4)+2]->SetTitle((elementStrNames[i] + " Integral Graph Mean").c_str());

        meanSeperateGraphs[(i*4)+3] = new TGraphErrors(cycles, timeArr, (cycleMeanResults[6]->getDoubleArrStorage())[i], zero, (cycleMeanResults[7]->getDoubleArrStorage())[i]);
        meanSeperateGraphs[(i*4)+3]->GetXaxis()->SetTitle("Time");
        meanSeperateGraphs[(i*4)+3]->GetYaxis()->SetTitle("Fit Value");
        meanSeperateGraphs[(i*4)+3]->SetTitle((elementStrNames[i] + " Single Integral Graph Mean").c_str());

    }
    delete [] zero;

    return meanSeperateGraphs;
}

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
}

//runs the cycles and puts the data of the integral method mean and the regular method mean into their respective arrays (dynamic)
resultStorage<Double_t>** Cycle::runSeperateMean()
{
    //setting inital values for time, 1= start const end moving, 2= end const start moving
    element->setTimeRunStart(x_start);
    element->setTimeRunEnd(x_stop);

    resultStorage<Double_t>** cycleSeperateMeanResult = new resultStorage<Double_t>* [8];
    cycleSeperateMeanResult[0] = new resultStorage<Double_t>(2, cycles, numElements); //Total regular
    cycleSeperateMeanResult[1] = new resultStorage<Double_t>(2, cycles, numElements); //Total regular error
    cycleSeperateMeanResult[2] = new resultStorage<Double_t>(2, cycles, numElements); //Total integral
    cycleSeperateMeanResult[3] = new resultStorage<Double_t>(2, cycles, numElements); //Total integral error
    cycleSeperateMeanResult[4] = new resultStorage<Double_t>(2, cycles, numElements); //Single regular
    cycleSeperateMeanResult[5] = new resultStorage<Double_t>(2, cycles, numElements); //Single regular error
    cycleSeperateMeanResult[6] = new resultStorage<Double_t>(2, cycles, numElements); //Single integral
    cycleSeperateMeanResult[7] = new resultStorage<Double_t>(2, cycles, numElements); //single integral error

    TH1D** multiRunHisto = decayChainRun->createRunResultHistos();
    TH1D** multiRunHistoSingle = decayChainRun->createRunResultHistosSingleElements();

    //running the cycle and putting histograms means in respective arrays
    for(int i = 0; i < cycles; i++)
    {
        decayChainRun->runNoChange();
        multiRunHisto = decayChainRun->fillRunResultHistos(multiRunHisto);
        multiRunHistoSingle = decayChainRun->fillRunResultHistosSingleElement(multiRunHistoSingle);
        for(int k = 0; k < numElements; k++)
        {
            (cycleSeperateMeanResult[0]->getDoubleArrStorage())[k][i] = multiRunHisto[(k*2)]->GetMean();
            //cout << "REG: " << multiRunHisto[(k*2)]->GetMean() << endl;
            (cycleSeperateMeanResult[1]->getDoubleArrStorage())[k][i] = multiRunHisto[(k*2)]->GetMeanError();
            //cout << "REG ERR: " << multiRunHisto[(k*2)]->GetMeanError() << endl << endl;
            (cycleSeperateMeanResult[2]->getDoubleArrStorage())[k][i] = multiRunHisto[(k*2)+1]->GetMean();
            //cout << "INTE: " << multiRunHisto[(k*2)+1]->GetMean() << endl;
            (cycleSeperateMeanResult[3]->getDoubleArrStorage())[k][i] = multiRunHisto[(k*2)+1]->GetMeanError();
            //cout << "INTE ERR: " << multiRunHisto[(k*2)+1]->GetMeanError() << endl << endl;
            (cycleSeperateMeanResult[4]->getDoubleArrStorage())[k][i] = multiRunHistoSingle[(k*2)]->GetMean();
            //cout << "REG: " << multiRunHistoSingle[(k*2)]->GetMean() << endl;
            (cycleSeperateMeanResult[5]->getDoubleArrStorage())[k][i] = multiRunHistoSingle[(k*2)]->GetMeanError();
            //cout << "REG ERR: " << multiRunHistoSingle[(k*2)]->GetMeanError() << endl << endl;
            (cycleSeperateMeanResult[6]->getDoubleArrStorage())[k][i] = multiRunHistoSingle[(k*2)+1]->GetMean();
            //cout << "INTE: " << multiRunHistoSingle[(k*2)+1]->GetMean() << endl;
            (cycleSeperateMeanResult[7]->getDoubleArrStorage())[k][i] = multiRunHistoSingle[(k*2)+1]->GetMeanError();
            //cout << "INTE ERR: " << multiRunHistoSingle[(k*2)+1]->GetMeanError() << endl << endl;
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
    return cycleSeperateMeanResult;
}

//does runs cycle for the single histogram (dynamic)
resultStorage<Double_t>** Cycle::runSeperateSingleGen()
{
    //setting inital values for time, 1 = start const end moving, 2 = start moving end const
    element->setTimeRunStart(x_start);
    element->setTimeRunEnd(x_stop);

    resultStorage<Double_t>** cycleSingleResult = new resultStorage<Double_t>* [8];
    cycleSingleResult[0] = new resultStorage<Double_t>(2, cycles, numElements); //Total regular
    cycleSingleResult[1] = new resultStorage<Double_t>(2, cycles, numElements); //Total regular error
    cycleSingleResult[2] = new resultStorage<Double_t>(2, cycles, numElements); //Total integral
    cycleSingleResult[3] = new resultStorage<Double_t>(2, cycles, numElements); //Total integral error
    cycleSingleResult[4] = new resultStorage<Double_t>(2, cycles, numElements); //Single regular
    cycleSingleResult[5] = new resultStorage<Double_t>(2, cycles, numElements); //Single regular error
    cycleSingleResult[6] = new resultStorage<Double_t>(2, cycles, numElements); //Single integral
    cycleSingleResult[7] = new resultStorage<Double_t>(2, cycles, numElements); //single integral error

    Double_t** values;
    //creates the histograms we will fit
    element->genRandomAlternate();
    element->genIntegralHisto();
    element->genIntegralSingleHistos();

    //running the cycle and puts the fit values in the arrays
    for(int i = 0; i < cycles; i++)
    {
        values = decayChainRun->runNoChange();
        for(int k = 0; k < numElements; k++)
        {
            (cycleSingleResult[0]->getDoubleArrStorage())[k][i] = values[2][k];
            //cout << "REG: " << values[2][k] << endl;
            (cycleSingleResult[1]->getDoubleArrStorage())[k][i] = values[0][k];
            //cout << "REG ERR: " << values[0][k] << endl << endl;
            (cycleSingleResult[2]->getDoubleArrStorage())[k][i] = values[3][k];
            //cout << "INTE: " << values[3][k] << endl;
            (cycleSingleResult[3]->getDoubleArrStorage())[k][i] = values[1][k];
            //cout << "INTE ERR: " << values[1][k] << endl << endl;
            (cycleSingleResult[4]->getDoubleArrStorage())[k][i] = values[8][k];
            //cout << "SINGLE REG: " << values[8][k] << endl;
            (cycleSingleResult[5]->getDoubleArrStorage())[k][i] = values[6][k];
            //cout << "SINGLE REG ERR: " << values[6][k] << endl << endl;
            (cycleSingleResult[6]->getDoubleArrStorage())[k][i] = values[9][k];
            //cout << "SINGLE INTE: " << values[9][k] << endl;
            (cycleSingleResult[7]->getDoubleArrStorage())[k][i] = values[7][k];
            //cout << "SINGLE INTE ERR: " << values[7][k] << endl << endl;
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

    return cycleSingleResult;
}
#endif