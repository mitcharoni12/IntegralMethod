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
    Int_t numRuns, numElements, inputHistoExecutionType, singleElementDataChoice;
    FitOption* fitOptions;                                  ///< Contains fit options for program.
    ElementFit* element;                                    ///< Element fit object. Used to do individual runs.
    ChainRunFitValues* batemanFitValues;                    ///< Stores fit values of different runs of the total Bateman histograms.
    ChainRunFitValues* integralFitValues;                   ///< Stores fit values of different runs of the total integral histograms.
    SingleChainRunFitValues* singleBatemanFitValues;        ///< Stores fit values of different runs of the single Bateman histograms.
    SingleChainRunFitValues* singleIntegralFitValues;       ///< Stores fit values of different runs of the single integral histograms.
    Double_t* eventsXAxis;                                  
    Double_t* timeFitEnd;                                   ///< Array containing fit end values for the input histogram.
    Double_t* runsXAxis;                                    ///< Array containing the runs indexes.
    Double_t *zero;                                         ///< Array of 0's, used for setting 0 error in the x values.
    string* elementNameStrs;                                ///< Contains element names for each element in the decay chain.
    TH1D** singleBatemanMultiRunFitHisto;                   ///< Histogram used to store the fit values for the multiple runs of the run class for the single bateman fit values.
    TH1D** singleIntegralMultiRunFitHisto;                  ///< Histogram used to store the fit values for the multiple runs of the run class for the single integral fit values
    TH1D** totalBatemanMultiRunFitHisto;                    ///< Histogram used to store the fit values for the multiple runs of the run class for the total bateman fit values
    TH1D** totalIntegralMultiRunFitHisto;                   ///< Histogram used to store the fit values for the multiple runs of the run class for the total integral fit values
    TGraphErrors** totalBatemanGraphs;                      ///< Stores graphs of fit values for the total Bateman histograms
    TGraphErrors** totalIntegralGraphs;                     ///< Stores graphs of fit values for the total integral histograms
    TGraphErrors** singleBatemanGraphs;                     ///< Stores graphs of fit values for the single Bateman histograms
    TGraphErrors** singleIntegralGraphs;                    ///< Stores graphs of fit values for the single integral histograms
    TCanvas* testCan;
public:
    Run(ElementFit* element);
    ~Run();
    //functions for single bateman data
    void CreateSingleBatemanMultiRunHistos();
    void DisplaySingleBatemanFitValuesGraphs(TCanvas** canvasArray);
    void DisplaySingleBatemanFitValuesHistos(TCanvas** canvasArray);
    void FillSingleBatemanMultiRunHistos();
    void GenSingleBatemanGraphsNoChange();
    void RunSingleBatemanRunsGenOnce(Int_t cycleIndex, Int_t runIndex);
    void RunSingleBatemanRunsNoChange(Int_t cycleIndex);
    //functions for single integral data
    void CreateSingleIntegralMultiRunHistos();
    void DisplaySingleIntegralFitValuesGraphs(TCanvas** canvasArray);
    void DisplaySingleIntegralFitValuesHistos(TCanvas** canvasArray);
    void FillSingleIntegralMultiRunHistos();
    void GenSingleIntegralGraphsNoChange();
    void RunSingleIntegralRunsGenOnce(Int_t cycleIndex, Int_t runIndex);
    void RunSingleIntegralRunsNoChange(Int_t cycleIndex);
    //functions for total bateman data
    void CreateTotalBatemanMultiRunHistos();
    void DisplayTotalBatemanFitValuesGraphs(TCanvas** canvasArray);
    void DisplayTotalBatemanFitValuesHistos(TCanvas** canvasArray);
    void FillTotalBatemanMultiRunHistos();
    void GenTotalBatemanGraphsNoChange();
    void RunTotalBatemanRunsGenOnce(Int_t cycleIndex, Int_t runIndex);
    void RunTotalBatemanRunsNoChange(Int_t cycleIndex);
    //functions for total integral data
    void CreateTotalIntegralMultiRunHistos();
    void DisplayTotalIntegralFitValuesGraphs(TCanvas** canvasArray);
    void DisplayTotalIntegralFitValuesHistos(TCanvas** canvasArray);
    void FillTotalIntegralMultiRunHistos();
    void GenTotalIntegralGraphsNoChange();
    void RunTotalIntegralRunsGenOnce(Int_t cycleIndex, Int_t runIndex);
    void RunTotalIntegralRunsNoChange(Int_t cycleIndex);

    //getter function
    ChainRunFitValues* GetBatemanFitValues(){return batemanFitValues;}
    ChainRunFitValues* GetIntegralFitValues(){return integralFitValues;}
    SingleChainRunFitValues* GetSingleBatemanFitValues(){return singleBatemanFitValues;}
    SingleChainRunFitValues* GetSingleIntegralFitValues(){return singleIntegralFitValues;}
    TH1D** GetSingleBatemanMultiRunFitHisto(){return singleBatemanMultiRunFitHisto;}
    TH1D** GetSingleIntegralMultiRunFitHisto(){return singleIntegralMultiRunFitHisto;}
    TH1D** GetTotalBatemanMultiRunFitHisto(){return totalBatemanMultiRunFitHisto;}
    TH1D** GetTotalIntegralMultiRunFitHisto(){return totalIntegralMultiRunFitHisto;}
    string* GetElementStringNames(){return elementNameStrs;}
    Int_t GetNumRuns(){return numRuns;}
    //setter function
    void setNumRuns(Int_t numRuns){this->numRuns = numRuns;}
};

/// \brief Constructor for Run.
///
/// Constructor for the run class. Dynamical allocation for the data storage objects for the fit values of the runs happens here.
Run::Run(ElementFit* element)
{
    //getting parameters for the program execution type
    this->fitOptions = element->GetFitOptions();
    this->numRuns = fitOptions->GetNumRuns();
    this->element = element;
    this->elementNameStrs = fitOptions->GetElementNames();
    this->timeFitEnd = fitOptions->GetTimeLengthArr();
    this->inputHistoExecutionType = fitOptions->GetInputHistoExecutionType();
    this->singleElementDataChoice = fitOptions->GetSingleElementDataChoice();
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
    //histos used to store fit values of the different runs
    totalBatemanMultiRunFitHisto = new TH1D* [numElements];
    totalIntegralMultiRunFitHisto = new TH1D* [numElements];
    if(singleElementDataChoice == 2)
    {
        singleBatemanMultiRunFitHisto = new TH1D* [numElements];
        singleIntegralMultiRunFitHisto = new TH1D* [numElements];
    }
    //used to store fit values for the different runs
    batemanFitValues = new ChainRunFitValues(numElements, numRuns);
    integralFitValues = new ChainRunFitValues(numElements, numRuns);
    if(singleElementDataChoice == 2)
    {
        singleBatemanFitValues = new SingleChainRunFitValues(numElements, numRuns);
        singleIntegralFitValues = new SingleChainRunFitValues(numElements, numRuns);
    }
    //used to plot fit values for different runs
    totalBatemanGraphs = new TGraphErrors* [numElements];
    totalIntegralGraphs = new TGraphErrors* [numElements];
    if(singleElementDataChoice == 2)
    {
        singleBatemanGraphs = new TGraphErrors* [numElements];
        singleIntegralGraphs = new TGraphErrors* [numElements];
    }
    //testCan = new TCanvas("testCanvas", "testCanvas", 500, 500);
}

Run::~Run()
{
    delete batemanFitValues;
    delete integralFitValues;
    delete [] runsXAxis;
    delete [] zero;
    delete [] totalBatemanGraphs;
    delete [] totalIntegralGraphs;
    delete [] totalBatemanMultiRunFitHisto;
    delete [] totalIntegralMultiRunFitHisto;
    if(singleElementDataChoice == 2)
    {
        delete singleBatemanFitValues;
        delete singleIntegralFitValues;
        delete [] singleBatemanGraphs;
        delete [] singleIntegralGraphs;
        delete [] singleBatemanMultiRunFitHisto;
        delete [] singleIntegralMultiRunFitHisto;
    }
}

/// \brief Dynamically allocates histograms for the storage of the fit values for the multiple runs of the single bateman histograms.
void Run::CreateSingleBatemanMultiRunHistos()
{
    Double_t parameterValue;

    //creates the histograms
    for(int i = 0; i < numElements; i++)
    {
        parameterValue = TMath::LogE()/(element->GetElementParameters(i));
        singleBatemanMultiRunFitHisto[i] = new TH1D((elementNameStrs[i] + " Fit Result Single Bateman Histo").c_str(), (elementNameStrs[i] + " Fit Result Single Bateman Histo").c_str(), 500, parameterValue*0, parameterValue*2.5);
    }
}

/// \brief Dynamically allocates histograms for the storage of the fit values for the multiple runs of the single integral histograms.
void Run::CreateSingleIntegralMultiRunHistos()
{
    Double_t parameterValue;

    //creates the histograms
    for(int i = 0; i < numElements; i++)
    {
        parameterValue = TMath::LogE()/(element->GetElementParameters(i));
        singleIntegralMultiRunFitHisto[i] = new TH1D((elementNameStrs[i] + " Fit Result Single Integral Histo").c_str(), (elementNameStrs[i] + " Fit Result Single Integral Histo").c_str(), 500, parameterValue*0, parameterValue*2.5);
    }
}

/// \brief Dynamically allocates histograms for the storage of the fit values for the multiple runs of the total bateman histograms.
void Run::CreateTotalBatemanMultiRunHistos()
{
    Double_t parameterValue;

    //creates the histograms
    for(int i = 0; i < numElements; i++)
    {
        parameterValue = TMath::LogE()/(element->GetElementParameters(i));
        //totalBatemanMultiRunFitHisto[i] = new TH1D((elementNameStrs[i] + " Fit Result Total Bateman Histo").c_str(), (elementNameStrs[i] + " Fit Result Total Bateman Histo").c_str(), 500, parameterValue*0, parameterValue*2.5);
        totalBatemanMultiRunFitHisto[i] = new TH1D(("Fitted Bateman Method Values").c_str(), (elementNameStrs[i] + " Fit Result Total Bateman Histo").c_str(), 500, parameterValue*0, parameterValue*2.5);
        totalBatemanMultiRunFitHisto[i]->GetXaxis()->SetTitle("Fitted Half-Life(S)");
        totalBatemanMultiRunFitHisto[i]->GetYaxis()->SetTitle("Counts");
    }
}

/// \brief Dynamically allocates histograms for the storage of the fit values for the multiple runs of the total integral histograms.
void Run::CreateTotalIntegralMultiRunHistos()
{
    Double_t parameterValue;

    //creates the histograms
    for(int i = 0; i < numElements; i++)
    {
        parameterValue = TMath::LogE()/(element->GetElementParameters(i));
        //totalIntegralMultiRunFitHisto[i] = new TH1D((elementNameStrs[i] + " Fit Result Total Integral Histo").c_str(), (elementNameStrs[i] + " Fit Result Total Integral Histo").c_str(), 500, parameterValue*0, parameterValue*2.5);
        totalIntegralMultiRunFitHisto[i] = new TH1D(("Fitted Integral Method Values").c_str(), (elementNameStrs[i] + " Fit Result Total Integral Histo").c_str(), 500, parameterValue*0, parameterValue*2.5);
        totalIntegralMultiRunFitHisto[i]->GetXaxis()->SetTitle("Fitted Half-Life(S)");
        totalIntegralMultiRunFitHisto[i]->GetYaxis()->SetTitle("Counts");
    }
}


/// \brief Displays a graph of the fit values of the multiple runs for the single bateman histograms.
void Run::DisplaySingleBatemanFitValuesGraphs(TCanvas** canvasArray)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(3);
        singleBatemanGraphs[i]->Draw();
    }
}

/// \brief Displays a graph of the fit values of the multiple runs for the single integral histograms.
void Run::DisplaySingleIntegralFitValuesGraphs(TCanvas** canvasArray)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(4);
        singleIntegralGraphs[i]->Draw();
    }
}

/// \brief Displays a graph of the fit values of the multiple runs for the total bateman histograms.
void Run::DisplayTotalBatemanFitValuesGraphs(TCanvas** canvasArray)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(1);
        totalBatemanGraphs[i]->Draw();
    }
}

/// \brief Displays a graph of the fit values of the multiple runs for the total integral histograms.
void Run::DisplayTotalIntegralFitValuesGraphs(TCanvas** canvasArray)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(2);
        totalIntegralGraphs[i]->Draw();
    }
}

/// \brief Displays a histogram of the fit values of the multiple runs for the single bateman histograms.
void Run::DisplaySingleBatemanFitValuesHistos(TCanvas** canvasArray)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(3);
        singleBatemanMultiRunFitHisto[i]->Draw();
    }
}

/// \brief Displays a histogram of the fit values of the multiple runs for the single integral histograms.
void Run::DisplaySingleIntegralFitValuesHistos(TCanvas** canvasArray)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(4);
        singleIntegralMultiRunFitHisto[i]->Draw();
    }
}

/// \brief Displays a histogram of the fit values of the multiple runs for the total bateman histograms.
void Run::DisplayTotalBatemanFitValuesHistos(TCanvas** canvasArray)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(1);
        totalBatemanMultiRunFitHisto[i]->Draw();
    }
}

/// \brief Displays a histogram of the fit values of the multiple runs for the total integral histograms.
void Run::DisplayTotalIntegralFitValuesHistos(TCanvas** canvasArray)
{
    for(int i = 0; i < numElements; i++)
    {
        canvasArray[i]->cd(2);
        totalIntegralMultiRunFitHisto[i]->Draw();
    }
}

/// \brief Fills the data storage histograms with fitted single bateman half lives from the runs.
void Run::FillSingleBatemanMultiRunHistos()
{
    for(int i = 0; i < numElements; i++)
    {
        singleBatemanMultiRunFitHisto[i]->Reset("ICES");
        for(int j = 0; j < numRuns; j++)
        {
            singleBatemanMultiRunFitHisto[i]->Fill(singleBatemanFitValues->GetAnHalfLife(j, i, i));
        }
    }
}

/// \brief Fills the data storage histograms with fitted single integral half lives from the runs.
void Run::FillSingleIntegralMultiRunHistos()
{
    for(int i = 0; i < numElements; i++)
    {
        singleIntegralMultiRunFitHisto[i]->Reset("ICES");
        for(int j = 0; j < numRuns; j++)
        {
            singleIntegralMultiRunFitHisto[i]->Fill(singleIntegralFitValues->GetAnHalfLife(j, i, i));
        }
    }
}

/// \brief Fills the data storage histograms with fitted total bateman half lives from the runs.
void Run::FillTotalBatemanMultiRunHistos()
{
    //delete later
    Int_t failedRuns = 0;
    for(int i = 0; i < numElements; i++)
    {
        totalBatemanMultiRunFitHisto[i]->Reset("ICES");
        for(int j = 0; j < numRuns; j++)
        {
            //IMPORTANT TO CHANGE LATER
            if(batemanFitValues->GetAnHalfLife(j, i) < -0.00001 || batemanFitValues->GetAnHalfLife(j, i) > 0.00001)
            {
                totalBatemanMultiRunFitHisto[i]->Fill(batemanFitValues->GetAnHalfLife(j, i));
            }else{
                failedRuns++;
            }
        }
        cout << "Element " << i << " Bateman mean: " << totalBatemanMultiRunFitHisto[i]->GetMean() << endl;
    }
        cout << "Bateman failed runs: " << failedRuns/numElements << endl;
}

/// \brief Fills the data storage histograms with fitted total integral half lives from the runs.
void Run::FillTotalIntegralMultiRunHistos()
{
    //delete later
    Int_t failedRuns = 0;
    for(int i = 0; i < numElements; i++)
    {
        totalIntegralMultiRunFitHisto[i]->Reset("ICES");
        for(int j = 0; j < numRuns; j++)
        {
            //IMPORTANT TO CHANGE LATER
            if(integralFitValues->GetAnHalfLife(j, i) < -0.00001 || integralFitValues->GetAnHalfLife(j, i) > 0.00001)
            {
                totalIntegralMultiRunFitHisto[i]->Fill(integralFitValues->GetAnHalfLife(j, i));
            }else{
                failedRuns++;
            }
        }
        cout << "Element " << i << " integral mean: " << totalIntegralMultiRunFitHisto[i]->GetMean() << endl;
    }
        cout << "Integral failed runs " << failedRuns/numElements << endl;
}

/// \brief Dynamically allocates the graphs to display the fit values of the single bateman histograms for the different runs.
void Run::GenSingleBatemanGraphsNoChange()
{
    for(int j = 0; j < numElements; j++)
    {
        singleBatemanGraphs[j]= new TGraphErrors(numRuns, runsXAxis, singleBatemanFitValues->GetHalfLifeArr(j, j), zero, singleBatemanFitValues->GetHalfLifeErrorArr(j, j));
        singleBatemanGraphs[j]->GetXaxis()->SetTitle("Runs");
        singleBatemanGraphs[j]->GetYaxis()->SetTitle("Bateman Fit(S)");
        singleBatemanGraphs[j]->SetTitle((elementNameStrs[j] + " bateman Fit Single Element(S)").c_str());
        singleBatemanGraphs[j]->SetName((elementNameStrs[j] + " bateman_Fit_Single_Element(S)").c_str());
    }
}

/// \brief Dynamically allocates the graphs to display the fit values of the single integral histograms for the different runs.
void Run::GenSingleIntegralGraphsNoChange()
{
    for(int j = 0; j < numElements; j++)
    {
        singleIntegralGraphs[j] = new TGraphErrors(numRuns, runsXAxis, singleIntegralFitValues->GetHalfLifeArr(j, j), zero, singleIntegralFitValues->GetHalfLifeErrorArr(j, j));
        singleIntegralGraphs[j]->GetXaxis()->SetTitle("Runs");
        singleIntegralGraphs[j]->GetYaxis()->SetTitle("Integral Fit Single Element(S)");
        singleIntegralGraphs[j]->SetTitle((elementNameStrs[j] + " Integral Fit Single Element(S)").c_str());
        singleIntegralGraphs[j]->SetName((elementNameStrs[j] + " Integral_Fit_Single_Element(S)").c_str());
    }
}

/// \brief Dynamically allocates the graphs to display the fit values of the total bateman histograms for the different runs.
void Run::GenTotalBatemanGraphsNoChange() 
{
    for(int j = 0; j < numElements; j++)
    {
        totalBatemanGraphs[j]= new TGraphErrors(numRuns, runsXAxis, batemanFitValues->GetHalfLifeArr(j), zero, batemanFitValues->GetHalfLifeErrorArr(j));
        totalBatemanGraphs[j]->GetXaxis()->SetTitle("Runs");
        totalBatemanGraphs[j]->GetYaxis()->SetTitle("Bateman Fit(S)");
        totalBatemanGraphs[j]->SetTitle((elementNameStrs[j] + " Bateman Fit(S)").c_str());
        totalBatemanGraphs[j]->SetName((elementNameStrs[j] + " bateman_Fit(S)").c_str());
    }
}

/// \brief Dynamically allocates the graphs to display the fit values of the total integral histograms for the different runs.
void Run::GenTotalIntegralGraphsNoChange() 
{
    for(int j = 0; j < numElements; j++)
    {
        totalIntegralGraphs[j] = new TGraphErrors(numRuns, runsXAxis, integralFitValues->GetHalfLifeArr(j), zero, integralFitValues->GetHalfLifeErrorArr(j));
        totalIntegralGraphs[j]->GetXaxis()->SetTitle("Runs");
        totalIntegralGraphs[j]->GetYaxis()->SetTitle("Integral Fit(S)");
        totalIntegralGraphs[j]->SetTitle((elementNameStrs[j] + " Integral Fit(S)").c_str());
        totalIntegralGraphs[j]->SetName((elementNameStrs[j] + " Integral_Fit(S)").c_str());
    }
}

/// \brief Preforms the multiple runs of the program for the single bateman data not changing anything setting between runs and only generating a single set of histograms to fit with for each cycle.
///
/// Only used with the cycle program execution type. It would not make sense to do this with the run program execution type becuase then it would be the same as just generating and fitting a single histogram.
void Run::RunSingleBatemanRunsGenOnce(Int_t cycleIndex, Int_t runIndex)
{
    SingleElementFitValues* singleTempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    element->FitSingleBatemanHistos(cycleIndex, 0);

    //get single bateman fit parameters
    singleTempFitParameters = element->GetSingleBatemanFitValues();

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
}

/// \brief Preforms the multiple runs of the program for the single integral data not changing anything setting between runs and only generating a single set of histograms to fit with for each cycle.
///
/// Only used with the cycle program execution type. It would not make sense to do this with the run program execution type becuase then it would be the same as just generating and fitting a single histogram.
void Run::RunSingleIntegralRunsGenOnce(Int_t cycleIndex, Int_t runIndex)
{
    SingleElementFitValues* singleTempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    element->FitSingleIntegralHistos(cycleIndex, 0);

    //get single integral fit parameters
    singleTempFitParameters = element->GetSingleIntegralFitValues();

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

/// \brief Preforms the multiple runs of the program for the total bateman data not changing anything setting between runs and only generating a single set of histograms to fit with for each cycle.
///
/// Only used with the cycle program execution type. It would not make sense to do this with the run program execution type becuase then it would be the same as just generating and fitting a single histogram.
void Run::RunTotalBatemanRunsGenOnce(Int_t cycleIndex, Int_t runIndex)
{
    ChainFitValues* tempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    element->FitTotalBatemanHistos(cycleIndex, 0);

    tempFitParameters = element->GetBatemanFitValues();

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
}

/// \brief Preforms the multiple runs of the program for the total integral data not changing anything setting between runs and only generating a single set of histograms to fit with for each cycle.
///
/// Only used with the cycle program execution type. It would not make sense to do this with the run program execution type becuase then it would be the same as just generating and fitting a single histogram.
void Run::RunTotalIntegralRunsGenOnce(Int_t cycleIndex, Int_t runIndex)
{
    ChainFitValues* tempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    element->FitTotalIntegralHistos(cycleIndex, 0);

    //gets total integral fit parameters
    tempFitParameters = element->GetIntegralFitValues();

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
}

/// \brief Preforms the multiple runs of the program for the single bateman data not changing anything setting between runs.
void Run::RunSingleBatemanRunsNoChange(Int_t cycleIndex)
{
    SingleElementFitValues* singleTempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    for(int j = 0; j < numRuns; j++)
    {
        //generates random data and fits it. Then extract the fit parametes
        element->FitSingleBatemanHistos(cycleIndex, j);

        //get single bateman fit parameters
        singleTempFitParameters = element->GetSingleBatemanFitValues();

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
    }
}

/// \brief Preforms the multiple runs of the program for the single integral data not changing anything setting between runs.
void Run::RunSingleIntegralRunsNoChange(Int_t cycleIndex)
{
    SingleElementFitValues* singleTempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    for(int j = 0; j < numRuns; j++)
    {
        //generates random data and fits it. Then extract the fit parametes
        element->FitSingleIntegralHistos(cycleIndex, j);

        //get single integral fit parameters
        singleTempFitParameters = element->GetSingleIntegralFitValues();

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

/// \brief Preforms the multiple runs of the program for the total bateman data not changing anything setting between runs.
void Run::RunTotalBatemanRunsNoChange(Int_t cycleIndex)
{
    ChainFitValues* tempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    for(int j = 0; j < numRuns; j++)
    {
        //generates random data and fits it. Then extract the fit parametes
        element->FitTotalBatemanHistos(cycleIndex, j);
        //gets total bateman fit parameters
        tempFitParameters = element->GetBatemanFitValues();
        
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
            cout << "Element: " << (i+1) << " run: " << (j+1) << " HalfLife: " << tempHalfLife << endl;
            batemanFitValues->SetAnHalfLifeError(j, i, tempHalfLifeError);
        }
        cout << endl;
    }
}

/// \brief Preforms the multiple runs of the program for the total integral data not changing anything setting between runs.
void Run::RunTotalIntegralRunsNoChange(Int_t cycleIndex)
{
    ChainFitValues* tempFitParameters;
    Double_t tempN0;
    Double_t tempN0Error;
    Double_t tempHalfLife;
    Double_t tempHalfLifeError;

    for(int j = 0; j < numRuns; j++)
    {
        //generates random data and fits it. Then extract the fit parametes
        element->FitTotalIntegralHistos(cycleIndex, j);
        //gets total integral fit parameters
        tempFitParameters = element->GetIntegralFitValues();

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
            cout << "Element: " << (i+1) << " run: " << (j+1) << " HalfLife: " << tempHalfLife << endl;
            integralFitValues->SetAnHalfLifeError(j, i, tempHalfLifeError);
        }
        cout << endl;
    }
}

#endif