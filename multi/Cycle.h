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
        FitOption* fitOptions;                                  ///< Contains fit options for program.
        Int_t numCycles, numElements;
        Int_t* binNumArr;                                       ///< Contains the bin number between cycles.
        Double_t* timeArr;                                      ///< Contains the time run end between cycles
        Run* decayChainRun;                                     ///< Run objects, used for doing the multiple runs of a cycle
        ElementFit* element;                                    ///< Element fit object. Used to do individual runs.
        string* elementStrNames;                                ///< Contains the names of the different elements in the decay chain.
        ChainRunFitValues* batemanFitValues;                    ///< Contains the average fit values for each individual cycle for the total bateman fit values. If its a single histogram generation type then it stores the different values fit between cycles.
        ChainRunFitValues* integralFitValues;                   ///< Contains the average fit values for each individual cycle for the total integral fit values. If its a single histogram generation type then it stores the different values fit between cycles.
        ChainRunFitValues* meanDifferenceValues;                ///< Contains the difference between the total integral and total bateman average fit values for each individual cycle for the total bateman fit values.(bateman - integral)
        SingleChainRunFitValues* singleBatemanFitValues;        ///< Contains the average fit values for each individual cycle for the single bateman fit values. If its a single histogram generation type then it stores the different values fit between cycles.
        SingleChainRunFitValues* singleIntegralFitValues;       ///< Contains the average fit values for each individual cycle for the single integral fit values. If its a single histogram generation type then it stores the different values fit between cycles.
        SingleChainRunFitValues* singleMeanDifferenceValues;    ///< Contains the difference between the single integral and single bateman average fit values for each individual cycle for the total bateman fit values.(bateman - integral)
        TGraphErrors** batemanFitMeanGraphs;                    ///< Graph contains the fit values between cycles for the total bateman histograms. If its a single histogram generation then its just the fits between the cycles of the single histogram.
        TGraphErrors** integralFitMeanGraphs;                   ///< Graph contains the fit values between cycles for the total integral histograms. If its a single histogram generation then its just the fits between the cycles of the single histogram.
        TGraphErrors** singleBatemanFitMeanGraphs;              ///< Graph contains the fit values between cycles for the single bateman histograms. If its a single histogram generation then its just the fits between the cycles of the single histogram.
        TGraphErrors** singleIntegralFitMeanGraphs;             ///< Graph contains the fit values between cycles for the single integral histograms. If its a single histogram generation then its just the fits between the cycles of the single histogram.
        TGraphErrors** fitDifferenceGraphs;                     ///< Graph contains the difference between the total bateman and total integral fit values. (bateman - integral)
        TGraphErrors** singleFitDifferenceGraphs;               ///< Graph contains the difference between the single bateman and single integral fit values. (bateman - integral)
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


/// \brief Constructor for Cycle class.
///
/// Constructor for Cycle class. Dynamic allocation of data storage objects for the class happen here.
Cycle::Cycle(Run* decayChainRun, ElementFit* element)
{
    this->fitOptions = element->getFitOptions();
    this->numCycles = fitOptions->GetNumCycles();
    this->element = element;
    this->decayChainRun = decayChainRun;
    numElements = fitOptions->GetNumElements();
    binNumArr = fitOptions->GetBinNumArr();
    elementStrNames = fitOptions->GetElementNames();
    batemanFitValues = new ChainRunFitValues(numElements, numCycles);
    integralFitValues = new ChainRunFitValues(numElements, numCycles);
    meanDifferenceValues = new ChainRunFitValues(numElements, numCycles);
    singleBatemanFitValues = new SingleChainRunFitValues(numElements, numCycles);
    singleIntegralFitValues = new SingleChainRunFitValues(numElements, numCycles);
    singleMeanDifferenceValues = new SingleChainRunFitValues(numElements, numCycles);
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

/// \brief Displays graphs for difference in mean of the bateman and integral fits.
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

/// \brief Displays the graphs for the mean of the bateman and integral fits.
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

/// \brief Generates graphs for rebinning between cycles and difference of Bateman and integral fits.
///
/// Dynamically allocates graphs for the rebinning between cycles and the difference of the Bateman and integral fits. After dynamic allocation the graphs are filled with the correct data.
void Cycle::genMeanDifferenceGraphsRebin()
{
    //create arrays for zero error and the bin numbers between each cycle
    Double_t* zero = new Double_t[numCycles];
    Double_t* binNumArrDoubles = new Double_t [numCycles];
    for(int i = 0; i < numCycles; i++)
    {
        zero[i] = 0.0f;
    }
    for(int i = 0; i < numCycles; i++)
    {
        binNumArrDoubles[i] = (Double_t) binNumArr[i];
    }

    //dynamically allocate and fill the graphs
    for(int i = 0; i < numElements; i++)
    {
        fitDifferenceGraphs[i] = new TGraphErrors(numCycles, binNumArrDoubles, meanDifferenceValues->GetHalfLifeArr(i), zero, meanDifferenceValues->GetHalfLifeErrorArr(i));
        fitDifferenceGraphs[i]->GetXaxis()->SetTitle("Number Bins");
        fitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit(s)");
        fitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Bateman Histo Mean Difference").c_str());

        singleFitDifferenceGraphs[i] = new TGraphErrors(numCycles, binNumArrDoubles, singleMeanDifferenceValues->GetHalfLifeArr(i, i), zero, singleMeanDifferenceValues->GetHalfLifeErrorArr(i, i));
        singleFitDifferenceGraphs[i]->GetXaxis()->SetTitle("Number Bins");
        singleFitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit(s)");
        singleFitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Single Histo Mean Difference").c_str());
    }
    delete [] zero;
    delete [] binNumArrDoubles;
}

/// \brief generates graphs for changing fit time and the difference between the Bateman and integral fits
///
/// Dynamically allocates and fills graphs for different fit times between runs for the difference between the Bateman and integral fits.
void Cycle::genMeanDifferenceGraphsTimeDifference()
{
    //dynamic allocation of zero error arrays
    Double_t* zero = new Double_t[numCycles];
    for(int i = 0; i < numCycles; i++)
    {
        zero[i] = 0.0f;
    }

    //dynamically allocate and fill the graphs
    for(int i = 0; i < numElements; i++)
    {
        fitDifferenceGraphs[i] = new TGraphErrors(numCycles, timeArr, meanDifferenceValues->GetHalfLifeArr(i), zero, meanDifferenceValues->GetHalfLifeErrorArr(i));
        fitDifferenceGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        fitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit(s)");
        fitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Bateman Histo Mean Difference").c_str());

        singleFitDifferenceGraphs[i] = new TGraphErrors(numCycles, timeArr, singleMeanDifferenceValues->GetHalfLifeArr(i, i), zero, singleMeanDifferenceValues->GetHalfLifeErrorArr(i, i));
        singleFitDifferenceGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        singleFitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit");
        singleFitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Single Histo Mean Difference(s)").c_str());
    }
    delete [] zero;
}

/// \brief Calculates the mean difference between Bateman and integral fits and puts them in their respective array.
void Cycle::genSingleMeanDifference()
{
    Double_t differenceHalfLife, differenceHalfLifeError;

    //calculates difference and error
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k < numCycles; k++)
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

/// \brief Dynamically allocates and fills the graphs for the mean fit values for the bateman and integral fits.
void Cycle::genSeperateMeanGraphsRebin()
{
    Double_t* tempFitVals, *tempFitErrors;
    Double_t* zero = new Double_t[numCycles];
    Double_t* binNumArrDoubles = new Double_t [numCycles];

    //dynamic allocation of zero error array and bin number between cycles array
    for(int i = 0; i < numCycles; i++)
    {
        zero[i] = 0.0f;
    }
    for(int i = 0; i < numCycles; i++)
    {
        binNumArrDoubles[i] = (Double_t) binNumArr[i];
    }

    //dynamically allocates graphs then fills them with values
    for(int i = 0; i < numElements; i++)
    {
        tempFitVals = batemanFitValues->GetHalfLifeArr(i);
        tempFitErrors = batemanFitValues->GetHalfLifeErrorArr(i);
        batemanFitMeanGraphs[i] = new TGraphErrors(numCycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        batemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        batemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        batemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Bateman Graph Mean").c_str());

        tempFitVals = singleBatemanFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleBatemanFitValues->GetHalfLifeErrorArr(i, i);
        singleBatemanFitMeanGraphs[i] = new TGraphErrors(numCycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        singleBatemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        singleBatemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        singleBatemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single bateman Graph Mean").c_str());

        tempFitVals = integralFitValues->GetHalfLifeArr(i);
        tempFitErrors = integralFitValues->GetHalfLifeErrorArr(i);
        integralFitMeanGraphs[i] = new TGraphErrors(numCycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        integralFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        integralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        integralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Integral Graph Mean").c_str());

        tempFitVals = singleIntegralFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleIntegralFitValues->GetHalfLifeErrorArr(i, i);
        singleIntegralFitMeanGraphs[i] = new TGraphErrors(numCycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        singleIntegralFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        singleIntegralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(S)");
        singleIntegralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Integral Graph Mean").c_str());
    }

    delete [] zero;
    delete [] binNumArrDoubles;
}

/// \brief Dynamically allocates and fills the graphs for different fit times between runs.
void Cycle::genSeperateMeanGraphsTimeDifference()
{
    Double_t* tempFitVals, *tempFitErrors;
    Double_t* zero = new Double_t[numCycles];

    //zero error array
    for(int i = 0; i < numCycles; i++)
    {
        zero[i] = 0.0f;
    }

    //dynamically allocates and fills graph
    for(int i = 0; i < numElements; i++)
    {
        tempFitVals = batemanFitValues->GetHalfLifeArr(i);
        tempFitErrors = batemanFitValues->GetHalfLifeErrorArr(i);
        batemanFitMeanGraphs[i] = new TGraphErrors(numCycles, timeArr, tempFitVals, zero, tempFitErrors);
        batemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        batemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        batemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Bateman Graph Mean").c_str());

        tempFitVals = singleBatemanFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleBatemanFitValues->GetHalfLifeErrorArr(i, i);
        singleBatemanFitMeanGraphs[i] = new TGraphErrors(numCycles, timeArr, tempFitVals, zero, tempFitErrors);
        singleBatemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        singleBatemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        singleBatemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Bateman Graph Mean").c_str());

        tempFitVals = integralFitValues->GetHalfLifeArr(i);
        tempFitErrors = integralFitValues->GetHalfLifeErrorArr(i);
        integralFitMeanGraphs[i] = new TGraphErrors(numCycles, timeArr, tempFitVals, zero, tempFitErrors);
        integralFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        integralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        integralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Integral Graph Mean").c_str());

        tempFitVals = singleIntegralFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleIntegralFitValues->GetHalfLifeErrorArr(i, i);
        singleIntegralFitMeanGraphs[i] = new TGraphErrors(numCycles, timeArr, tempFitVals, zero, tempFitErrors);
        singleIntegralFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        singleIntegralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        singleIntegralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Integral Graph Mean").c_str());
    }

    delete [] zero;
}

/// \brief Runs the Cycles and changes the number of bins between cycles. Bateman and integral fits subtracted and that value is stored.
///
/// Runs the Cycles and changes the bins between cycles. Once a single cycle is run the difference between the Bateman and integral fit is taken and that is stored in the corresponding data object.
void Cycle::runDifferenceMeanRebin()
{
    Double_t differenceHalfLife, differenceHalfLifeError;

    TH1D** multiRunHisto = decayChainRun->createRunResultHistos();
    TH1D** multiRunSingleHistos = decayChainRun->createRunResultHistosSingleElements();
    
    //elements in array stored in order, bateman, single
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->runNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs and returns that here.
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

/// \brief Runs the Cycles and changes the time fit between cycles. Bateman and integral fits subtracted and that value is stored.
///
/// Runs the Cycles and changes the time fit between cycles. Once a single cycle is run the difference between the Bateman and integral fit is taken and that is stored in the corresponding data object.
void Cycle::runDifferenceMeanTimeDifference()
{
    Double_t differenceHalfLife, differenceHalfLifeError;

    TH1D** multiRunHisto = decayChainRun->createRunResultHistos();
    TH1D** multiRunSingleHistos = decayChainRun->createRunResultHistosSingleElements();
    
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->runNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs and returns that here.
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

/// \brief Runs the Cycles and changes the number of bins between cycles.
///
/// Runs the Cycles and changes the bins between cycles. Once a single cycle is run the value is taken and stored in the corresponding data object.
void Cycle::runSeperateMeanRebin()
{
    Double_t tempMeanVal;
    Double_t tempErrorVal;

    TH1D** multiRunHisto = decayChainRun->createRunResultHistos();
    TH1D** multiRunHistoSingle = decayChainRun->createRunResultHistosSingleElements();

    //running the cycle and putting histograms means in respective arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->runNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs and returns that here.
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

/// \brief Runs the Cycles and changes the time fit between cycles.
///
/// Runs the Cycles and changes the time fit between cycles. Once a single cycle is run the value is taken and stored in the corresponding data object.
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
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->runNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs and returns that here.
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
        cout << elementStrNames[i] << " AVERAGE BATEMAN TOTAL: " << runningTotalBateman[i] / numCycles << endl;
        cout << elementStrNames[i] << " AVERAGE INTEGRAL TOTAL: " << runningTotalIntegral[i] / numCycles << endl;
        cout << elementStrNames[i] << " AVERAGE BATEMAN SINGLE: " << runningSingleBateman[i] / numCycles << endl;
        cout << elementStrNames[i] << " AVERAGE INTEGRAL SINGLE: " << runningSingleIntegral[i] / numCycles << endl << endl;
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
/// \brief Runs the cycles for the generation and fitting of a single histogram.
///
/// Runs the cycles for the generation of a single histogram, changing whatever variable(time fit) between each cycle.
void Cycle::runSeperateSingleGen()
{
    Double_t N0;
    Double_t N0Error;
    Double_t halfLife;
    Double_t halfLifeError;

    ChainRunFitValues* tempVals;
    SingleChainRunFitValues* singleTempVals;

    //running the cycle and puts the fit values in the arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
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