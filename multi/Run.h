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

/// Class will handle multiple runs of the program. Used extensively in the Cycle class.
class Run{
private:
    Int_t numRuns, numElements, eventDecrement;
    FitOption* fitOptions;                                  ///< Contains fit options for program.
    ElementFit* element;                                    ///< Element fit object. Used to do individual runs.
    ChainRunFitValues* batemanFitValues;                    ///< Stores fit values of different runs of the total Bateman histograms.
    ChainRunFitValues* integralFitValues;                   ///< Stores fit values of different runs of the total integral histograms.
    SingleChainRunFitValues* singleBatemanFitValues;        ///< Stores fit values of different runs of the single Bateman histograms.
    SingleChainRunFitValues* singleIntegralFitValues;       ///< Stores fit values of different runs of the single integral histograms.
    Double_t* eventsXAxis;                                  
    Double_t* runsXAxis;                                    ///< Array containing the runs indexes.
    Double_t *zero;                                         ///< Array of 0's, used for setting 0 error in the x values.
    string* elementNameStrs;                                ///< Contains element names for each element in the decay chain.
    TH1D** multiRunResultHistograms;                        ///< Histograms used to store the fit values
    TGraphErrors** totalBatemanGraphs;                      ///< Stores graphs of fit values for the total Bateman histograms
    TGraphErrors** totalIntegralGraphs;                     ///< Stores graphs of fit values for the total integral histograms
    TGraphErrors** singleBatemanGraphs;                     ///< Stores graphs of fit values for the single Bateman histograms
    TGraphErrors** singleIntegralGraphs;                    ///< Stores graphs of fit values for the single integral histograms
    TCanvas* testCan;
    //helper functions
    Double_t getMaxElement(Double_t* arr);
    Double_t getMinElement(Double_t* arr);
    void genIntegralMeanGraphs();
    void genBatemanMeanGraphs();
public:
    Run(ElementFit* element);
    ~Run();
    void runNoChangeGenOnce(Int_t cycleIndex, Int_t runIndex);
    //void genGraphsEventChange();
    void genGraphsNoChange();
    void genGraphsNoChangeSingleElement();
    TH1D** createRunResultHistos();
    TH1D** createRunResultHistosSingleElements();
    TH1D** fillRunResultHistos(TH1D** multiRunResultHistograms);
    TH1D** fillRunResultHistosSingleElement(TH1D** multiRunResultHistogramsSingleElement);
    void displayMultiRunResultGraphs(TCanvas** canvasArray);
    void displayMultiRunResultHistos(TCanvas** canvasArray, TH1D** multiRunResultHistograms);
    void runNoChange(Int_t cycleIndex);
    //getter function
    ChainRunFitValues* getBatemanFitValues(){return batemanFitValues;}
    ChainRunFitValues* getIntegralFitValues(){return integralFitValues;}
    SingleChainRunFitValues* getSingleBatemanFitValues(){return singleBatemanFitValues;}
    SingleChainRunFitValues* getSingleIntegralFitValues(){return singleIntegralFitValues;}
    TH1D** getMultiRunResultHistos(){return multiRunResultHistograms;}
    string* getElementStringNames(){return elementNameStrs;}
    Int_t getNumRuns(){return numRuns;}
    //setter function
    void setNumRuns(Int_t numRuns){this->numRuns = numRuns;}
};

/// \brief Constructor for Run.
///
/// Constructor for the run class. Dynamical allocation for the data storage objects for the fit values of the runs happens here.
Run::Run(ElementFit* element)
{
    //getting parameters for the program execution type
    this->fitOptions = element->getFitOptions();
    this->numRuns = fitOptions->GetNumRuns();
    this->eventDecrement = fitOptions->GetEventDecrement();
    this->element = element;
    this->elementNameStrs = fitOptions->GetElementNames();
    numElements = fitOptions->GetNumElements();

    //dynamical allocation for arrays containing events of general
    runsXAxis = new Double_t [numRuns];
    for(int i = 0; i < numRuns; i++)
    {
        runsXAxis[i] = i+1;
    }
    zero = new Double_t [numRuns];
    for(int i = 0; i < numRuns; i++)
    {
        zero[i] = 0.0f;
    }
    for(int i = 0; i < numRuns; i++)
    {
        runsXAxis[i] = i+1;
    }
    //used to store fit values for the different runs
    batemanFitValues = new ChainRunFitValues(numElements, numRuns);
    integralFitValues = new ChainRunFitValues(numElements, numRuns);
    singleBatemanFitValues = new SingleChainRunFitValues(numElements, numRuns);
    singleIntegralFitValues = new SingleChainRunFitValues(numElements, numRuns);
    //used to plot fit values for different runs
    totalBatemanGraphs = new TGraphErrors* [numElements];
    totalIntegralGraphs = new TGraphErrors* [numElements];
    singleBatemanGraphs = new TGraphErrors* [numElements];
    singleIntegralGraphs = new TGraphErrors* [numElements];
    //testCan = new TCanvas("testCanvas", "testCanvas", 500, 500);
}

Run::~Run()
{
    delete batemanFitValues;
    delete integralFitValues;
    delete singleBatemanFitValues;
    delete singleIntegralFitValues;
    delete [] runsXAxis;
    delete [] totalBatemanGraphs;
    delete [] totalIntegralGraphs;
    delete [] singleBatemanGraphs;
    delete [] singleIntegralGraphs;
}

/// \brief Dynamically allocates histograms for the storage of the fit values
///
/// Dynamically allocates histograms to store fit values of total histograms for multiple runs of the program.
/// Main use if for getting average fit values for x amount of runs.
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

/// \brief Dynamically allocates histograms for the storage of the fit values
///
/// Dynamically allocates histograms to store fit values of single histograms for multiple runs of the program.
/// Main use if for getting average fit values for x amount of runs.
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

/// \brief Displays graphs for the fit values of the multiple runs.
///
/// Displays graphs for the fit values of multiple runs. Will display each fit values individually in a graph like format.
/// Useful for seeing if fitting is conistent. Graphs displayed grouped by element. Displays data for the 4 different fit types.
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

/// \brief Displays histogram for the fit values of multiple runs.
///
/// Displays the histograms for the fit values of multile runs. Displays the fit values for each element in a histogram.
/// Useful for seeing the average fitted value. Histogram displayed grouped by element. Display data for 4 different fit types.
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

/// \brief Fills the data storage histograms with fitted data from the runs.
///
/// Fills the data storage histograms with fitted data from the runs. Only retreives half lives values and retreives the fit data for the total histograms.
TH1D** Run::fillRunResultHistos(TH1D** multiRunResultHistograms)
{
    for(int i = 0; i < numElements; i++)
    {
        multiRunResultHistograms[(i*2)]->Reset("ICES");
        multiRunResultHistograms[(i*2)+1]->Reset("ICES");
        for(int j = 0; j < numRuns; j++)
        {
            multiRunResultHistograms[(i*2)]->Fill(batemanFitValues->GetAnHalfLife(j, i));
            multiRunResultHistograms[(i*2)+1]->Fill(integralFitValues->GetAnHalfLife(j, i));
        }
    }

    return multiRunResultHistograms;
}

/// \brief Fills the data storage histograms with fitted data from the runs.
///
/// Fills the data storage histograms with fitted data from the runs. Only retreives half lives values and retreives the fit data for the single histograms.
TH1D** Run::fillRunResultHistosSingleElement(TH1D** multiRunResultHistogramsSingleElement)
{
    for(int i = 0; i < numElements; i++)
    {
        multiRunResultHistogramsSingleElement[(i*2)]->Reset("ICES");
        multiRunResultHistogramsSingleElement[(i*2)+1]->Reset("ICES");
        for(int j = 0; j < numRuns; j++)
        {
            multiRunResultHistogramsSingleElement[(i*2)]->Fill(singleBatemanFitValues->GetAnHalfLife(j, i, i));
            multiRunResultHistogramsSingleElement[(i*2)+1]->Fill(singleIntegralFitValues->GetAnHalfLife(j, i, i));
        }
    }

    return multiRunResultHistogramsSingleElement;
}

/*
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
*/

/// \brief Dynamically allocates the graphs to display the fit values of the total histograms for the different runs.
///
/// Dynamically allocates the graphs to display the fit values of the total histograms for the different runs. Retreives the array of fit values and creates the graphs.
/// Then do things like set graph name. Only do this for half lives. Graphs are generated in the conideration nothing is chaning between the runs.
void Run::genGraphsNoChange() 
{
    for(int j = 0; j < numElements; j++)
    {
        totalBatemanGraphs[j]= new TGraphErrors(numRuns, runsXAxis, batemanFitValues->GetHalfLifeArr(j), zero, batemanFitValues->GetHalfLifeErrorArr(j));
        totalBatemanGraphs[j]->GetXaxis()->SetTitle("Runs");
        totalBatemanGraphs[j]->GetYaxis()->SetTitle("Bateman Fit(S)");
        totalBatemanGraphs[j]->SetTitle((elementNameStrs[j] + " Bateman Fit(S)").c_str());
        totalBatemanGraphs[j]->SetName((elementNameStrs[j] + " bateman_Fit(S)").c_str());

        totalIntegralGraphs[j] = new TGraphErrors(numRuns, runsXAxis, integralFitValues->GetHalfLifeArr(j), zero, integralFitValues->GetHalfLifeErrorArr(j));
        totalIntegralGraphs[j]->GetXaxis()->SetTitle("Runs");
        totalIntegralGraphs[j]->GetYaxis()->SetTitle("Integral Fit(S)");
        totalIntegralGraphs[j]->SetTitle((elementNameStrs[j] + " Integral Fit(S)").c_str());
        totalIntegralGraphs[j]->SetName((elementNameStrs[j] + " Integral_Fit(S)").c_str());
    }
}

/// \brief Dynamically allocates the graphs to display the fit values of the single histograms for the different runs.
///
/// Dynamically allocates the graphs to display the fit values of the single histograms for the different runs. Retreives the array of fit values and creates the graphs.
/// Then do things like set graph name. Only do this for half lives. Graphs are generated in the conideration nothing is chaning between the runs.
void Run::genGraphsNoChangeSingleElement()
{
    for(int j = 0; j < numElements; j++)
    {
        singleBatemanGraphs[j]= new TGraphErrors(numRuns, runsXAxis, singleBatemanFitValues->GetHalfLifeArr(j, j), zero, singleBatemanFitValues->GetHalfLifeErrorArr(j, j));
        singleBatemanGraphs[j]->GetXaxis()->SetTitle("Runs");
        singleBatemanGraphs[j]->GetYaxis()->SetTitle("Bateman Fit(S)");
        singleBatemanGraphs[j]->SetTitle((elementNameStrs[j] + " bateman Fit Single Element(S)").c_str());
        singleBatemanGraphs[j]->SetName((elementNameStrs[j] + " bateman_Fit_Single_Element(S)").c_str());

        singleIntegralGraphs[j] = new TGraphErrors(numRuns, runsXAxis, singleIntegralFitValues->GetHalfLifeArr(j, j), zero, singleIntegralFitValues->GetHalfLifeErrorArr(j, j));
        singleIntegralGraphs[j]->GetXaxis()->SetTitle("Runs");
        singleIntegralGraphs[j]->GetYaxis()->SetTitle("Integral Fit Single Element(S)");
        singleIntegralGraphs[j]->SetTitle((elementNameStrs[j] + " Integral Fit Single Element(S)").c_str());
        singleIntegralGraphs[j]->SetName((elementNameStrs[j] + " Integral_Fit_Single_Element(S)").c_str());
    }
}

/// \brief Gets max value out of array
Double_t Run::getMaxElement(Double_t* arr)
{
    Double_t hold = arr[0];
    for(int i = 1; i < numRuns; i++)
    {
        if(hold < arr[i])
        {
            hold  = arr[i];
        }
    }
    return hold;
}

/// \brief Gets minimum value out of array
Double_t Run::getMinElement(Double_t* arr)
{
    Double_t hold = arr[0];
    for(int i = 1; i < numRuns; i++)
    {
        if(hold > arr[i])
        {
            hold  = arr[i];
        }
    }
    return hold;
}

/// \brief Preforms the multiple runs of the program not changing anything setting between runs.
///
/// Preforms the multiple runs of the program not changing anything setting between runs. The function first fits a histogram at a specific run and cycle index.
/// It then stores that data in the repsective data objects.
void Run::runNoChange(Int_t cycleIndex)
{
    //dynamic array
    ChainFitValues* tempFitParameters;
    SingleElementFitValues* singleTempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    for(int j = 0; j < numRuns; j++)
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


/// \brief Preforms the multiple runs of the program not changing anything setting between runs and only generating a single set of histograms to fit with for each cycle.
///
/// Preforms the multiple runs of the program not changing anything setting between runs and only generating a single set of histograms to fit with for each cycle.
/// Only used with the cycle program execution type. It would not make sense to do this with the run program execution type becuase then it would be the same as just generating and fitting a single histogram.
void Run::runNoChangeGenOnce(Int_t cycleIndex, Int_t runIndex)
{
    ChainFitValues* tempFitParameters;
    SingleElementFitValues* singleTempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    element->fitHistos(cycleIndex, runIndex);

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