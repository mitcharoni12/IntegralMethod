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
    Int_t numCycles, numElements, inputHistoExecutionType, singleElementDataChoice;
    bool displayFitAverages;                                ///< Determines if the fit averages for the different fit types will be displayed
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
    //functions for single element bateman and integral fits difference
    void DisplaySingleMeanDifferenceGraphs(TCanvas** canvasArr);
    void GenSingleMeanDifferenceGraphsRebin();
    void GenSingleMeanDifferenceGraphsTimeDifference();
    void GenSingleMeanDifference();
    void RunSingleDifferenceMeanCycles();
    //functions for total element bateman and integral fits difference
    void DisplayTotalMeanDifferenceGraphs(TCanvas** canvasArr);
    void GenTotalMeanDifferenceGraphsRebin();
    void GenTotalMeanDifferenceGraphsTimeDifference();
    void GenTotalMeanDifference();
    void RunTotalDifferenceMeanCycles();
    //functions for single bateman fits
    void DisplaySingleBatemanMeanGraphs(TCanvas** canvasArr);
    void GenSingleBatemanMeanGraphsRebin();
    void GenSingleBatemanMeanGraphsTimeDifference();
    void RunSingleBatemanCycles();
    void RunSingleBatemanCyclesSingleGen();
    //functions for single integral fits
    void DisplaySingleIntegralMeanGraphs(TCanvas** canvasArr);
    void GenSingleIntegralMeanGraphsRebin();
    void GenSingleIntegralMeanGraphsTimeDifference();
    void RunSingleIntegralCycles();
    void RunSingleIntegralCyclesSingleGen();
    //functions for total bateman fits
    void DisplayTotalBatemanMeanGraphs(TCanvas** canvasArr);
    void GenTotalBatemanMeanGraphsRebin();
    void GenTotalBatemanMeanGraphsTimeDifference();
    void RunTotalBatemanCycles();
    void RunTotalBatemanCyclesSingleGen();
    //functions for total integral fits
    void DisplayTotalIntegralMeanGraphs(TCanvas** canvasArr);
    void GenTotalIntegralMeanGraphsRebin();
    void GenTotalIntegralMeanGraphsTimeDifference();
    void RunTotalIntegralCycles();
    void RunTotalIntegralCyclesSingleGen();
};


/// \brief Constructor for Cycle class.
///
/// Constructor for Cycle class. Dynamic allocation of data storage objects for the class happen here.
Cycle::Cycle(Run* decayChainRun, ElementFit* element)
{
    this->fitOptions = element->GetFitOptions();
    this->numCycles = fitOptions->GetNumCycles();
    this->element = element;
    this->decayChainRun = decayChainRun;
    this->inputHistoExecutionType = fitOptions->GetInputHistoExecutionType();
    this->singleElementDataChoice = fitOptions->GetSingleElementDataChoice();
    numElements = fitOptions->GetNumElements();
    binNumArr = fitOptions->GetBinNumArr();
    displayFitAverages = fitOptions->GetDisplayFitAverages();
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

/// \brief Displays graphs for difference in mean of the signle bateman and single integral fits.
void Cycle::DisplaySingleMeanDifferenceGraphs(TCanvas** canvasArr)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(2);
        singleFitDifferenceGraphs[i]->Draw();
    }
}

/// \brief Displays graphs for difference in mean of the total bateman and total integral fits.
void Cycle::DisplayTotalMeanDifferenceGraphs(TCanvas** canvasArr)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(1);
        fitDifferenceGraphs[i]->Draw();
    }
}

/// \brief Displays the graphs for the mean of the single bateman and fits.
void Cycle::DisplaySingleBatemanMeanGraphs(TCanvas** canvasArr)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(3);
        singleBatemanFitMeanGraphs[i]->Draw();
        canvasArr[i]->cd(4);
        singleIntegralFitMeanGraphs[i]->Draw();
    }
}

/// \brief Displays the graphs for the mean of the single integral and fits.
void Cycle::DisplaySingleIntegralMeanGraphs(TCanvas** canvasArr)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(4);
        singleIntegralFitMeanGraphs[i]->Draw();
    }
}

/// \brief Displays the graphs for the mean of the total bateman and fits.
void Cycle::DisplayTotalBatemanMeanGraphs(TCanvas** canvasArr)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(1);
        batemanFitMeanGraphs[i]->Draw();
    }
}

/// \brief Displays the graphs for the mean of the total integral and fits.
void Cycle::DisplayTotalIntegralMeanGraphs(TCanvas** canvasArr)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArr[i]->cd(2);
        integralFitMeanGraphs[i]->Draw();
    }
}

/// \brief Generates graphs for rebinning between cycles and difference of single bateman and single integral fits(bateman - integral).
void Cycle::GenSingleMeanDifferenceGraphsRebin()
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
        singleFitDifferenceGraphs[i] = new TGraphErrors(numCycles, binNumArrDoubles, singleMeanDifferenceValues->GetHalfLifeArr(i, i), zero, singleMeanDifferenceValues->GetHalfLifeErrorArr(i, i));
        singleFitDifferenceGraphs[i]->GetXaxis()->SetTitle("Number Bins");
        singleFitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit(s)");
        singleFitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Single Histo Mean Difference").c_str());
    }
    delete [] zero;
    delete [] binNumArrDoubles;
}

/// \brief Generates graphs for rebinning between cycles and difference of total bateman and total integral fits(bateman - integral).
void Cycle::GenTotalMeanDifferenceGraphsRebin()
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
    }
    delete [] zero;
    delete [] binNumArrDoubles;
}

/// \brief generates graphs for changing fit times between the cycles for the difference between the single bateman and single integral fits
void Cycle::GenSingleMeanDifferenceGraphsTimeDifference()
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
        singleFitDifferenceGraphs[i] = new TGraphErrors(numCycles, timeArr, singleMeanDifferenceValues->GetHalfLifeArr(i, i), zero, singleMeanDifferenceValues->GetHalfLifeErrorArr(i, i));
        singleFitDifferenceGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        singleFitDifferenceGraphs[i]->GetYaxis()->SetTitle("Bateman Fit - Integral Fit");
        singleFitDifferenceGraphs[i]->SetTitle((elementStrNames[i] + " Single Histo Mean Difference(s)").c_str());
    }
    delete [] zero;
}

/// \brief generates graphs for changing fit times between the cycles for the difference between the total Bateman and total integral fits
void Cycle::GenTotalMeanDifferenceGraphsTimeDifference()
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
    }
    delete [] zero;
}

/// \brief Calculates the mean difference between the single bateman and single integral fits and puts them in their respective array.
void Cycle::GenSingleMeanDifference()
{
    Double_t differenceHalfLife, differenceHalfLifeError;

    //calculates difference and error
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k < numCycles; k++)
        {
            differenceHalfLife = singleBatemanFitValues->GetAnHalfLife(k, i, i) - singleIntegralFitValues->GetAnHalfLife(k, i, i);
            singleMeanDifferenceValues->SetAnHalfLife(k, i, i, differenceHalfLife);

            differenceHalfLifeError = sqrt(TMath::Power(singleBatemanFitValues->GetAnHalfLifeError(k, i, i),2) + TMath::Power(singleIntegralFitValues->GetAnHalfLifeError(k, i, i),2));
            singleMeanDifferenceValues->SetAnHalfLifeError(k, i, i, differenceHalfLifeError);
        }
    }
}

/// \brief Calculates the mean difference between the total bateman and total integral fits and puts them in their respective array.
void Cycle::GenTotalMeanDifference()
{
    Double_t differenceHalfLife, differenceHalfLifeError;

    //calculates difference and error
    for(int i = 0; i < numElements; i++)
    {
        for(int k = 0; k < numCycles; k++)
        {
            differenceHalfLife = batemanFitValues->GetAnHalfLife(k, i) - integralFitValues->GetAnHalfLife(k, i);
            meanDifferenceValues->SetAnHalfLife(k, i, differenceHalfLife);

            differenceHalfLifeError = sqrt(TMath::Power(batemanFitValues->GetAnHalfLifeError(k, i),2) + TMath::Power(integralFitValues->GetAnHalfLifeError(k, i),2));
            meanDifferenceValues->SetAnHalfLifeError(k, i, differenceHalfLifeError);
        }
    }
}

/// \brief Dynamically allocates and fills the graphs for the mean fit values for the single bateman fits.
void Cycle::GenSingleBatemanMeanGraphsRebin()
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
        tempFitVals = singleBatemanFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleBatemanFitValues->GetHalfLifeErrorArr(i, i);
        singleBatemanFitMeanGraphs[i] = new TGraphErrors(numCycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        singleBatemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        singleBatemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        singleBatemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single bateman Graph Mean").c_str());
    }

    delete [] zero;
    delete [] binNumArrDoubles;
}

/// \brief Dynamically allocates and fills the graphs for the mean fit values for the single integral fits.
void Cycle::GenSingleIntegralMeanGraphsRebin()
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

/// \brief Dynamically allocates and fills the graphs for the mean fit values for the total bateman fits.
void Cycle::GenTotalBatemanMeanGraphsRebin()
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
    }

    delete [] zero;
    delete [] binNumArrDoubles;
}

/// \brief Dynamically allocates and fills the graphs for the mean fit values for the total integral fits.
void Cycle::GenTotalIntegralMeanGraphsRebin()
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
        tempFitVals = integralFitValues->GetHalfLifeArr(i);
        tempFitErrors = integralFitValues->GetHalfLifeErrorArr(i);
        integralFitMeanGraphs[i] = new TGraphErrors(numCycles, binNumArrDoubles, tempFitVals, zero, tempFitErrors);
        integralFitMeanGraphs[i]->GetXaxis()->SetTitle("Bin Number");
        integralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        integralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Integral Graph Mean").c_str());
    }

    delete [] zero;
    delete [] binNumArrDoubles;
}

/// \brief Dynamically allocates and fills the graphs for different fit times between cycles for the single bateman fit.
void Cycle::GenSingleBatemanMeanGraphsTimeDifference()
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
        tempFitVals = singleBatemanFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleBatemanFitValues->GetHalfLifeErrorArr(i, i);
        singleBatemanFitMeanGraphs[i] = new TGraphErrors(numCycles, timeArr, tempFitVals, zero, tempFitErrors);
        singleBatemanFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        singleBatemanFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        singleBatemanFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Bateman Graph Mean").c_str());
    }

    delete [] zero;
}

/// \brief Dynamically allocates and fills the graphs for different fit times between cycles for the single integral fit.
void Cycle::GenSingleIntegralMeanGraphsTimeDifference()
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
        tempFitVals = singleIntegralFitValues->GetHalfLifeArr(i, i);
        tempFitErrors = singleIntegralFitValues->GetHalfLifeErrorArr(i, i);
        singleIntegralFitMeanGraphs[i] = new TGraphErrors(numCycles, timeArr, tempFitVals, zero, tempFitErrors);
        singleIntegralFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        singleIntegralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        singleIntegralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Single Integral Graph Mean").c_str());
    }

    delete [] zero;
}

/// \brief Dynamically allocates and fills the graphs for different fit times between cycles for the single bateman fit.
void Cycle::GenTotalBatemanMeanGraphsTimeDifference()
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
    }

    delete [] zero;
}

/// \brief Dynamically allocates and fills the graphs for different fit times between cycles for the total integral fit.
void Cycle::GenTotalIntegralMeanGraphsTimeDifference()
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
        tempFitVals = integralFitValues->GetHalfLifeArr(i);
        tempFitErrors = integralFitValues->GetHalfLifeErrorArr(i);
        integralFitMeanGraphs[i] = new TGraphErrors(numCycles, timeArr, tempFitVals, zero, tempFitErrors);
        integralFitMeanGraphs[i]->GetXaxis()->SetTitle("Time(s)");
        integralFitMeanGraphs[i]->GetYaxis()->SetTitle("Fit Value(s)");
        integralFitMeanGraphs[i]->SetTitle((elementStrNames[i] + " Integral Graph Mean").c_str());
    }

    delete [] zero;
}

/// \brief Runs the Cycles and changes the number of bins between cycles. single bateman and single integral fits subtracted and that value is stored(bateman - integral).
void Cycle::RunSingleDifferenceMeanCycles()
{
    Double_t differenceHalfLife, differenceHalfLifeError;
    TH1D** tempBatemanHistos, **tempIntegralHistos;

    decayChainRun->CreateSingleBatemanMultiRunHistos();
    decayChainRun->CreateSingleIntegralMultiRunHistos();
    
    //elements in array stored in order, bateman, single
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunSingleBatemanRunsNoChange(i);
        decayChainRun->RunSingleIntegralRunsNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs and returns that here.
        decayChainRun->FillSingleBatemanMultiRunHistos();
        decayChainRun->FillSingleIntegralMultiRunHistos();

        tempBatemanHistos = decayChainRun->GetSingleBatemanMultiRunFitHisto();
        tempIntegralHistos = decayChainRun->GetSingleIntegralMultiRunFitHisto();

        //does mean difference calculation for each element in decay chain
        for(int k = 0; k < numElements; k++)
        {
            differenceHalfLife = tempBatemanHistos[k]->GetMean() - tempIntegralHistos[k]->GetMean();
            singleMeanDifferenceValues->SetAnHalfLife(i, k, k, differenceHalfLife);

            differenceHalfLifeError = sqrt(pow(tempBatemanHistos[k]->GetMeanError(),2) + pow(tempIntegralHistos[k]->GetMeanError(),2));
            singleMeanDifferenceValues->SetAnHalfLifeError(i, k, k, differenceHalfLifeError);
        }
    }
}

/// \brief Runs the Cycles and changes the number of bins between cycles. total bateman and total integral fits subtracted and that value is stored(bateman - integral).
void Cycle::RunTotalDifferenceMeanCycles()
{
    Double_t differenceHalfLife, differenceHalfLifeError;
    TH1D** tempBatemanHistos, **tempIntegralHistos;

    decayChainRun->CreateTotalBatemanMultiRunHistos();
    decayChainRun->CreateTotalIntegralMultiRunHistos();
    
    //elements in array stored in order, bateman, single
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunTotalBatemanRunsNoChange(i);
        decayChainRun->RunTotalIntegralRunsNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs and returns that here.
        decayChainRun->FillTotalBatemanMultiRunHistos();
        decayChainRun->FillTotalIntegralMultiRunHistos();

        tempBatemanHistos = decayChainRun->GetTotalBatemanMultiRunFitHisto();
        tempIntegralHistos = decayChainRun->GetTotalIntegralMultiRunFitHisto();

        //does mean difference calculation for each element in decay chain
        for(int k = 0; k < numElements; k++)
        {
            differenceHalfLife = tempBatemanHistos[k]->GetMean() - tempIntegralHistos[k]->GetMean();
            meanDifferenceValues->SetAnHalfLife(i, k, differenceHalfLife);

            differenceHalfLifeError = sqrt(pow(tempBatemanHistos[k]->GetMeanError(),2) + pow(tempIntegralHistos[k]->GetMeanError(),2));
            meanDifferenceValues->SetAnHalfLifeError(i, k, differenceHalfLifeError);
        }
    }
}

/// \brief Runs the multiple cycles of the program for multi source option for the single bateman fits.
void Cycle::RunSingleBatemanCycles()
{
    Double_t tempMeanVal, tempErrorVal;
    Double_t* runningSingleBateman;
    TH1D** tempFitHisto; 
    if(displayFitAverages)
    {
        runningSingleBateman = new Double_t [numElements];
        for(int i = 0; i < numElements; i++)
        {
            runningSingleBateman[i] = 0.0f;
        }
    }

    decayChainRun->CreateSingleBatemanMultiRunHistos();

    //running the cycle and putting histograms means in respective arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunSingleBatemanRunsNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs.
        decayChainRun->FillSingleBatemanMultiRunHistos();
        //gets the filled full of fit values
        tempFitHisto = decayChainRun->GetSingleBatemanMultiRunFitHisto();

        //moves histo means into respective data storages
        for(int k = 0; k < numElements; k++)
        {
            tempMeanVal = tempFitHisto[k]->GetMean();
            singleBatemanFitValues->SetAnHalfLife(i, k, k, tempMeanVal);
            singleBatemanFitValues->SetAnHalfLifeError(i, k, k, tempErrorVal);
            if(displayFitAverages)
            {
                runningSingleBateman[k] += tempMeanVal;
            }
            tempErrorVal = tempFitHisto[k]->GetMeanError();
        }
    }

    if(displayFitAverages)
    {
        for(int i = 0; i < numElements; i++)
        {
            cout << elementStrNames[i] << " AVERAGE BATEMAN SINGLE: " << runningSingleBateman[i] / numCycles << endl << endl;
        }

        delete [] runningSingleBateman;
    }
}

/// \brief Runs the multiple cycles of the program for multi source option for the single integral fits.
void Cycle::RunSingleIntegralCycles()
{
    Double_t tempMeanVal, tempErrorVal;
    Double_t* runningSingleIntegral;
    TH1D** tempFitHisto; 
    if(displayFitAverages)
    {
        runningSingleIntegral = new Double_t [numElements];
        for(int i = 0; i < numElements; i++)
        {
            runningSingleIntegral[i] = 0.0f;
        }
    }

    decayChainRun->CreateSingleIntegralMultiRunHistos();

    //running the cycle and putting histograms means in respective arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunSingleIntegralRunsNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs and returns that here.
        decayChainRun->FillSingleIntegralMultiRunHistos();
        //gets the filled full of fit values
        tempFitHisto = decayChainRun->GetSingleIntegralMultiRunFitHisto();

        //moves histo means into respective data storages
        for(int k = 0; k < numElements; k++)
        {
            tempMeanVal = tempFitHisto[k]->GetMean();
            singleIntegralFitValues->SetAnHalfLife(i, k, k, tempMeanVal);
            tempErrorVal = tempFitHisto[k]->GetMeanError();
            singleIntegralFitValues->SetAnHalfLifeError(i, k, k, tempErrorVal);
            if(displayFitAverages)
            {
                runningSingleIntegral[k] += tempMeanVal;
            }
        }
    }

    if(displayFitAverages)
    {
        for(int i = 0; i < numElements; i++)
        {
            cout << elementStrNames[i] << " AVERAGE INTEGRAL SINGLE: " << runningSingleIntegral[i] / numCycles << endl << endl;
        }
        delete [] runningSingleIntegral;
    }
}

/// \brief Runs the multiple cycles of the program for multi source option for the total bateman fits.
void Cycle::RunTotalBatemanCycles()
{
    Double_t tempMeanVal, tempErrorVal;
    Double_t* runningTotalBateman;
    TH1D** tempFitHisto; 
    if(displayFitAverages)
    {
        runningTotalBateman = new Double_t [numElements];
        for(int i = 0; i < numElements; i++)
        {
            runningTotalBateman[i] = 0.0f;
        }
    }

    decayChainRun->CreateTotalBatemanMultiRunHistos();

    //running the cycle and putting histograms means in respective arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunTotalBatemanRunsNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs and returns that here.
        decayChainRun->FillTotalBatemanMultiRunHistos();
        //gets the filled full of fit values
        tempFitHisto = decayChainRun->GetTotalBatemanMultiRunFitHisto();

        //moves histo means into respective data storages
        for(int k = 0; k < numElements; k++)
        {
            tempMeanVal = tempFitHisto[k]->GetMean();
            batemanFitValues->SetAnHalfLife(i, k, tempMeanVal);
            tempErrorVal = tempFitHisto[k]->GetMeanError();
            batemanFitValues->SetAnHalfLifeError(i, k, tempErrorVal);
            if(displayFitAverages)
            {
                runningTotalBateman[k] += tempMeanVal;
            }
        }
    }

    if(displayFitAverages)
    {
        for(int i = 0; i < numElements; i++)
        {
            cout << elementStrNames[i] << " AVERAGE BATEMAN TOTAL: " << runningTotalBateman[i] / numCycles << endl << endl;
        }
        delete [] runningTotalBateman;
    }
}

/// \brief Runs the multiple cycles of the program for multi source option for the total integral fits.
void Cycle::RunTotalIntegralCycles()
{
    Double_t tempMeanVal, tempErrorVal;
    Double_t *runningTotalIntegral;
    TH1D** tempFitHisto; 
    if(displayFitAverages)
    {
        runningTotalIntegral = new Double_t [numElements];
        for(int i = 0; i < numElements; i++)
        {
            runningTotalIntegral[i] = 0.0f;
        }
    }

    decayChainRun->CreateTotalIntegralMultiRunHistos();

    //running the cycle and putting histograms means in respective arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunTotalIntegralRunsNoChange(i);
        //fills a histogram with the fit values of the (numRuns) runs and returns that here.
        decayChainRun->FillTotalIntegralMultiRunHistos();
        //gets the filled full of fit values
        tempFitHisto = decayChainRun->GetTotalIntegralMultiRunFitHisto();

        //moves histo means into respective data storages
        for(int k = 0; k < numElements; k++)
        {
            tempMeanVal = tempFitHisto[k]->GetMean();
            integralFitValues->SetAnHalfLife(i, k, tempMeanVal);
            tempErrorVal = tempFitHisto[k]->GetMeanError();
            integralFitValues->SetAnHalfLifeError(i, k, tempErrorVal);
            if(displayFitAverages)
            {
                runningTotalIntegral[k] += tempMeanVal;
            }
        }
    }

    if(displayFitAverages)
    {
        for(int i = 0; i < numElements; i++)
        {
            cout << elementStrNames[i] << " AVERAGE INTEGRAL TOTAL: " << runningTotalIntegral[i] / numCycles << endl << endl;
        }
        delete [] runningTotalIntegral;
    }
}

/// \brief Runs the cycles for the generation and fitting of a single histogram for the single bateman fits.
void Cycle::RunSingleBatemanCyclesSingleGen()
{
    Double_t N0, N0Error, halfLife, halfLifeError;

    SingleChainRunFitValues* singleTempVals;

    //running the cycle and puts the fit values in the arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunSingleBatemanRunsGenOnce(i, 0);
        
        //gets single bateman fit data from fit
        singleTempVals = decayChainRun->GetSingleBatemanFitValues();

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
    }
}

/// \brief Runs the cycles for the generation and fitting of a single histogram for the single integral fits.
void Cycle::RunSingleIntegralCyclesSingleGen()
{
    Double_t N0, N0Error, halfLife, halfLifeError;

    SingleChainRunFitValues* singleTempVals;

    //running the cycle and puts the fit values in the arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunSingleIntegralRunsGenOnce(i, 0);
        
        //gets single bateman fit data from fit
        singleTempVals = decayChainRun->GetSingleIntegralFitValues();

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

/// \brief Runs the cycles for the generation and fitting of a single histogram for the total bateman fits.
void Cycle::RunTotalBatemanCyclesSingleGen()
{
    Double_t N0, N0Error, halfLife, halfLifeError;

    ChainRunFitValues* tempVals;

    //running the cycle and puts the fit values in the arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunTotalBatemanRunsGenOnce(i, 0);
        
        //gets total bateman fit data from fit
        tempVals = decayChainRun->GetBatemanFitValues();

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
    }
}

/// \brief Runs the cycles for the generation and fitting of a single histogram for the total integral fits.
void Cycle::RunTotalIntegralCyclesSingleGen()
{
    Double_t N0, N0Error, halfLife, halfLifeError;

    ChainRunFitValues* tempVals;

    //running the cycle and puts the fit values in the arrays
    for(int i = 0; i < numCycles; i++)
    {
        //runs (numRuns) runs
        decayChainRun->RunTotalIntegralRunsGenOnce(i, 0);
        
        //gets integral fit data from fit
        tempVals = decayChainRun->GetIntegralFitValues();

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
    }
}

#endif