#include <fstream>

#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TFile.h"

#include "ElementFit.h"
#include "Run.h"
#include "Cycle.h"
#include "ParameterValue.h"
#include "CycleCanvasHolder.h"
#include "SingleCycleCanvasHolder.h"
#include "FitFunction.h"

using namespace std;

Double_t IntegralDecaybyActivity(Double_t *x, Double_t *par);
Double_t BatemanDecaybyActivity(Double_t *x, Double_t *par);
void CreateFitFunctions(string* elementNames);

typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);
decayFunction* batemanFitFunctions;
decayFunction* integralFitFunctions;
Int_t numElements;

void fitMultiple();

int main()
{
    fitMultiple();
    return 0;
}

void fitMultiple()
{
    Int_t numBins, rebinBinInc, createFitFunctionsChoice, events, programExecutionType, histoSimulateChoice, rangeValueChoice,
          displayIndividualFitsChoice, rebinChoice, inputHistoExecutionType, inputHistoBinNum, singleElementDataChoice;
    Double_t timeRunEndSimulated, timeRunEndInput, timeRunStartInput, inputHistoTimeEnd, valueHolder, leaveOutStartBinsSim, leaveOutEndBinsSim, binWidth, leaveOutStartBinsInput, leaveOutEndBinsInput;
    FitOption* fitOptions = new FitOption();
    ElementFit* element;
    ifstream inFile;
    string rootFilePath, histogramName, elementDataFile;
    TFile* inputRootFile;
    TH1D* inputHistogram;

    //fitFunctions array passed in and used in element object
    /*
    decayFunction* fitFunctions = new decayFunction[numFitFunction];
    fitFunctions[0] = CSDecaybyActivity;
    fitFunctions[1] = CSDecaybyActivityIntegral;
    fitFunctions[2] = BADecaybyActivity;
    fitFunctions[3] = BADecaybyActivityIntegral;
    fitFunctions[4] = LADecaybyActivity;
    fitFunctions[5] = LADecaybyActivityIntegral;
    
    //arrays used to create the full decay function in a modular fassion in main
    batemanFitFunctions = new decayFunction [numElements];
    batemanFitFunctions[0] = CSDecaybyActivity;
    batemanFitFunctions[1] = BADecaybyActivity;
    batemanFitFunctions[2] = LADecaybyActivity;
    integralFitFunctions = new decayFunction [numElements];
    integralFitFunctions[0] = CSDecaybyActivityIntegral;
    integralFitFunctions[1] = BADecaybyActivityIntegral;
    integralFitFunctions[2] = LADecaybyActivityIntegral;
    */

    //READ IN DATA FROM SIMULATE.TXT
    inFile.open("simulate.txt");
    //gets numElements, events, bins, and timeRun
    inFile.ignore(256,':');
    inFile >> numElements;
    inFile.ignore(256,':');
    inFile >> events;
    inFile.ignore(256,':');
    inFile >> numBins;
    inFile.ignore(256,':');
    inFile >> inputHistoBinNum;
    inFile.ignore(256,':');
    inFile >> binWidth;
    inFile.ignore(256,':');
    inFile >> timeRunEndSimulated;
    inFile.ignore(256,':');
    inFile >> timeRunStartInput;
    inFile.ignore(256,':');
    inFile >> timeRunEndInput;
    inFile.ignore(256,':');
    inFile >> inputHistoTimeEnd;
    inFile.ignore(256,':');
    inFile >> singleElementDataChoice;
    //gets the fit generation choice
    inFile.ignore(256,':');
    inFile >> programExecutionType;
    inFile.ignore(256,':');
    inFile >> inputHistoExecutionType;
    inFile.ignore(256,':');
    inFile >> leaveOutStartBinsSim;
    inFile.ignore(256,':');
    inFile >> leaveOutEndBinsSim;
    inFile.ignore(256,':');
    inFile >> leaveOutStartBinsInput;
    inFile.ignore(256,':');
    inFile >> leaveOutEndBinsInput;
    inFile.ignore(256,':');
    inFile >> elementDataFile;
    inFile.ignore(256,':');
    inFile >> rootFilePath;
    inFile.ignore(256,':');
    inFile >> histogramName;
    inFile.close();

    //setting data read in from simulate.txt to the fitOption object.
    fitOptions->SetNumEvents(events);
    fitOptions->SetNumBins(numBins);
    fitOptions->SetBinWidth(binWidth);
    fitOptions->SetTimeRunEndSimulated(timeRunEndSimulated);
    fitOptions->SetTimeRunEndInput(timeRunEndInput);
    fitOptions->SetTimeRunStartInput(timeRunStartInput);
    fitOptions->SetProgramExecutionType(programExecutionType);
    fitOptions->SetLeaveOutStartBinsSim(leaveOutStartBinsSim);
    fitOptions->SetLeaveOutEndBinsSim(leaveOutEndBinsSim);
    fitOptions->SetInputHistoBinNum(inputHistoBinNum);
    fitOptions->SetLeaveOutStartBinsInput(leaveOutStartBinsInput);
    fitOptions->SetLeaveOutEndBinsInput(leaveOutEndBinsInput);
    fitOptions->SetNumElements(numElements);
    fitOptions->SetInputHistoExecutionType(inputHistoExecutionType);
    fitOptions->SetInputHistoTimeEnd(inputHistoTimeEnd);
    fitOptions->SetSingleElementDataChoice(singleElementDataChoice);
    Int_t numFitFunction = numElements*2;

    //open first option file
    inFile.open(elementDataFile.c_str());
    //create storage for the read in parameter values
    ParameterValue** paraVals = new ParameterValue* [numElements];
    //ENTERING VALUES FROM PARAMETER VALUES
    for(int i = 0; i < numElements; i++)
    {
        paraVals[i] = new ParameterValue();
    }
    string* elementNames = new string [numElements];
    for(int i = 0; i < numElements; i++)
    {
        //name of element
        inFile.ignore(256,':');
        inFile >> elementNames[i];
        
        //Gets Bateman N0 value and range
        inFile.ignore(256,':');
        inFile >> valueHolder;
        paraVals[i]->SetBatemanN0(valueHolder);
        inFile.ignore(256,':');
        inFile >> valueHolder;
        paraVals[i]->SetLowerRangeBatemanN0(valueHolder);
        inFile.ignore(256,':');
        inFile >> valueHolder;
        paraVals[i]->SetUpperRangeBatemanN0(valueHolder);
        inFile.ignore(256,':');
        inFile >> valueHolder;
        if(valueHolder == 1)
        {
            paraVals[i]->SetFixBatemanN0(true);
        }else
        {
            paraVals[i]->SetFixBatemanN0(false);
        }

        //Gets Bateman N0 value and range
        inFile.ignore(256,':');
        inFile >> valueHolder;
        paraVals[i]->SetIntegralN0(valueHolder);
        inFile.ignore(256,':');
        inFile >> valueHolder;
        paraVals[i]->SetLowerRangeIntegralN0(valueHolder);
        inFile.ignore(256,':');
        inFile >> valueHolder;
        paraVals[i]->SetUpperRangeIntegralN0(valueHolder);
        inFile.ignore(256,':');
        inFile >> valueHolder;
        if(valueHolder == 1)
        {
            paraVals[i]->SetFixIntegralN0(true);
        }else
        {
            paraVals[i]->SetFixIntegralN0(false);
        }

        //Gets Half life value and range
        inFile.ignore(256,':');
        inFile >> valueHolder;
        paraVals[i]->SetHalfLife(valueHolder);
        inFile.ignore(256,':');
        inFile >> valueHolder;
        paraVals[i]->SetLowerRangeHalfLife(valueHolder);
        inFile.ignore(256,':');
        inFile >> valueHolder;
        paraVals[i]->SetUpperRangeHalfLife(valueHolder);
        inFile.ignore(256,':');
        inFile >> valueHolder;
        if(valueHolder == 1)
        {
            paraVals[i]->SetFixHalfLife(true);
        }else
        {
            paraVals[i]->SetFixHalfLife(false);
        }
    }
    //gets option to create the fit functions if need be
    inFile.ignore(256,':');
    inFile >> createFitFunctionsChoice;
    fitOptions->SetElementNames(elementNames);
    inFile.close();

    //creates the fit functions in FitFunction.h
    if(createFitFunctionsChoice == 1)
    {
        CreateFitFunctions(elementNames);
        delete fitOptions;
        for(int i = 0; i < numElements; i++)
        {
            delete paraVals[i];
        }
        delete [] paraVals;
        delete [] elementNames;
        return;
    }
    FitFunction* fitFunctions = new FitFunction(numElements);
    batemanFitFunctions = fitFunctions->GetBatemanFitFunctions();
    integralFitFunctions = fitFunctions->GetIntegralFitFunctions();

    //calculates the decay constant values
    //note: lower range half life will produce upper range decay constant becuase of the formula for conversion
    for(int i = 0; i < numElements; i++)
    {
        valueHolder = (log(2)/paraVals[i]->GetLowerRangeHalfLife());
        paraVals[i]->SetUpperRangeDecayConst(valueHolder);
        valueHolder = (log(2)/paraVals[i]->GetUpperRangeHalfLife());
        paraVals[i]->SetLowerRangeDecayConst(valueHolder);
        valueHolder = (log(2)/paraVals[i]->GetHalfLife());
        paraVals[i]->SetDecayConst(valueHolder);
    }
    
    //switch statement for dealing with an input histogram
    switch(inputHistoExecutionType)
    {   //Case for pure histogram simulation.
        case 1:
        {
            break;
        }
        //Case for inputting histogram to do a Monte Carlo analysis
        case 2:
        {
            TH1D* inputHistogram;
            FitOption* inputHistoMonteFitOptions = fitOptions;
            ChainFitValues* integralFitValues;
            inputRootFile = new TFile(rootFilePath.c_str(), "READ");
            Double_t timeFitEnd, tempHalfLife;

            //if program cannot open root file
            if(inputRootFile->IsZombie())
            {
                cout << "Error reading file." << endl;
                delete fitOptions;
                for(int i = 0; i < numElements; i++)
                {
                    delete paraVals[i];
                }
                delete [] paraVals;
                delete [] elementNames;
                delete fitFunctions;
                return;
            }
            inputRootFile->GetObject(histogramName.c_str(), inputHistogram);

            inputHistoMonteFitOptions->SetNumRuns(1);
            inputHistoMonteFitOptions->SetNumCycles(1);
            element = new ElementFit(BatemanDecaybyActivity, IntegralDecaybyActivity, batemanFitFunctions, integralFitFunctions, paraVals, inputHistoMonteFitOptions, inputHistogram);

            timeFitEnd = inputHistoMonteFitOptions->GetTimeLengthArr()[0];
            element->FitTotalBatemanHistos(0, 0);
            element->FitTotalIntegralHistos(0, 0);
            element->DisplayTotalBatemanParameters();
            element->DisplayTotalIntegralParameters();
            integralFitValues = element->GetIntegralFitValues();

            for(int i = 0; i < numElements; i++)
            {
                tempHalfLife = integralFitValues->GetAnHalfLife(i);
                paraVals[i]->SetHalfLife(tempHalfLife);
            }

            delete element;
            fitOptions->SetSingleElementDataChoice(1);
            fitOptions->SetProgramExecutionType(2);
            cout << "\n\n\nSIMULATING DATA" << endl; 
            break;
        }
        //Case for changing fit time on input histogram.
        case 3:
        {
            Int_t numCycles, timeShiftType, displayIndividualFitsChoice, lowerCycleHistoIndex, upperCycleHistoIndex;
            Double_t timeInc;
            TH1D* inputHistogram;
            FitOption* inputHistoFitOptions = new FitOption();
            inputRootFile = new TFile(rootFilePath.c_str(), "READ");

            //if program cannot open root file
            if(inputRootFile->IsZombie())
            {
                cout << "Error reading file." << endl;
                delete fitOptions;
                for(int i = 0; i < numElements; i++)
                {
                    delete paraVals[i];
                }
                delete [] paraVals;
                delete [] elementNames;
                delete fitFunctions;
                return;
            }
            inputRootFile->GetObject(histogramName.c_str(), inputHistogram);

            inFile.open("inputHistoTimeChange.txt");
            inFile.ignore(256,':');
            inFile >> numCycles;
            inFile.ignore(256,':');
            inFile >> timeShiftType;
            inFile.ignore(256,':');
            inFile >> timeInc;
            inFile.ignore(256,':');
            inFile >> displayIndividualFitsChoice;
            inFile.ignore(256,':');
            inFile >> lowerCycleHistoIndex;
            inFile.ignore(256,':');
            inFile >> upperCycleHistoIndex;
            inFile.close();

            inputHistoFitOptions->SetNumRuns(1);
            inputHistoFitOptions->SetNumCycles(numCycles);
            inputHistoFitOptions->SetTimeShiftType(timeShiftType);
            inputHistoFitOptions->SetInputTimeInc(timeInc);
            inputHistoFitOptions->SetTimeRunEndInput(timeRunEndInput);
            inputHistoFitOptions->SetTimeRunStartInput(timeRunStartInput);
            inputHistoFitOptions->SetInputHistoBinNum(inputHistoBinNum);
            inputHistoFitOptions->SetLeaveOutStartBinsInput(leaveOutStartBinsInput);
            inputHistoFitOptions->SetLeaveOutEndBinsInput(leaveOutEndBinsInput);
            inputHistoFitOptions->SetNumElements(numElements);
            inputHistoFitOptions->SetInputHistoExecutionType(inputHistoExecutionType);
            inputHistoFitOptions->SetInputHistoTimeEnd(inputHistoTimeEnd);
            inputHistoFitOptions->SetElementNames(elementNames);
            inputHistoFitOptions->SetSingleElementDataChoice(1);

            element = new ElementFit(BatemanDecaybyActivity, IntegralDecaybyActivity, batemanFitFunctions, integralFitFunctions, paraVals, inputHistoFitOptions, inputHistogram);
            Run* run = new Run(element);
            Cycle* cycle = new Cycle(run, element);

            cycle->RunTotalBatemanCyclesSingleGen();
            cycle->RunTotalIntegralCyclesSingleGen();
            cycle->GenTotalBatemanMeanGraphsTimeDifference();
            cycle->GenTotalIntegralMeanGraphsTimeDifference();
            TCanvas** canvasArr = new TCanvas* [numElements];
            for(int i = 0; i < numElements; i++)
            {
                canvasArr[i] = new TCanvas((elementNames[i] + " Single Source Seperate Mean").c_str(), (elementNames[i] + " Single Source Seperate Mean").c_str(), 1100, 1100);
                canvasArr[i]->Divide(2,1,.02,.02);
            }
            cycle->DisplayTotalBatemanMeanGraphs(canvasArr);
            cycle->DisplayTotalIntegralMeanGraphs(canvasArr);
            delete [] canvasArr;

            if(displayIndividualFitsChoice == 1)
            {
                CycleCanvasHolder* batemanInputHistoCanvas = new CycleCanvasHolder(0, 0, lowerCycleHistoIndex, upperCycleHistoIndex, "Bateman Input Histogram");
                CycleCanvasHolder* integralInputHistoCanvas = new CycleCanvasHolder(0, 0, lowerCycleHistoIndex, upperCycleHistoIndex, "Integral Input Histogram");

                element->DrawTotalBatemanIndividualHistos(batemanInputHistoCanvas, 0, 0, lowerCycleHistoIndex, upperCycleHistoIndex);
                element->DrawTotalIntegralIndividualHistos(integralInputHistoCanvas, 0, 0, lowerCycleHistoIndex, upperCycleHistoIndex);

                delete batemanInputHistoCanvas;
                delete integralInputHistoCanvas;
            }

            delete inputHistoFitOptions;
            delete element;
            delete run;
            delete cycle;
            return;
        }
    }
    //switching program to simulate data now
    fitOptions->SetInputHistoExecutionType(1);

    //SWITCH FOR PROGRAM EXECUTION TYPE
    switch(fitOptions->GetProgramExecutionType())
    {
        //single run of histogam
        case 1:
        {
            fitOptions->SetMultiSourceChoice(true);
            element = new ElementFit(BatemanDecaybyActivity, IntegralDecaybyActivity, batemanFitFunctions, integralFitFunctions, paraVals, fitOptions);

            //fit total bateman and integral data
            element->FitTotalBatemanHistos(0, 0);
            element->FitTotalIntegralHistos(0, 0);

            //display total bateman and integral fitted parameters
            element->DisplayTotalBatemanParameters();
            element->DisplayTotalIntegralParameters();

            TCanvas* batemanTotalCanvas = new TCanvas("Total Bateman Histo", "Total Bateman Histo", 500, 500);
            TCanvas* integralTotalCanvas = new TCanvas("Total Integral Histo", "Total Integral Histo", 500, 500);

            //display total bateman and integral histograms
            element->DisplayTotalBatemanHisto(batemanTotalCanvas);
            element->DisplayTotalIntegralHisto(integralTotalCanvas);

            if(singleElementDataChoice == 2)
            {
                //fit single bateman and integral data
                element->FitSingleBatemanHistos(0, 0);
                element->FitSingleIntegralHistos(0, 0);

                //display single bateman and integral fitted parameters
                element->DisplaySingleBatemanParameters();
                element->DisplaySingleIntegralParameters();

                TCanvas** singleBatemanCanvases = new TCanvas* [numElements];
                TCanvas** singleIntegralCanvases = new TCanvas* [numElements];
                for(int i = 0; i < numElements; i++)
                {
                    singleBatemanCanvases[i] = new TCanvas((elementNames[i] + " Single Bateman Histo").c_str(), (elementNames[i] + " Single Bateman Histo").c_str(), 500, 500);
                    singleIntegralCanvases[i] = new TCanvas((elementNames[i] + " Single Integral Histo").c_str(), (elementNames[i] + " Single Integral Histo").c_str(), 500, 500);
                }

                //display total bateman and integral histogram
                element->DisplaySingleBatemanHistos(singleBatemanCanvases);
                element->DisplaySingleIntegralHisto(singleIntegralCanvases);

                delete [] singleBatemanCanvases;
                delete [] singleIntegralCanvases;
            }
            delete element;
        break;
        }
        //multi run of histogram
        case 2:
        {
            Int_t numRuns;
            int graphOrHistoChoice, lowerRunHistoIndex, upperRunHistoIndex;

            inFile.open("simulatedSingleCycle.txt");
            inFile.ignore(256,':');
            inFile >> numRuns;
            inFile.ignore(256,':');
            inFile >> graphOrHistoChoice;
            inFile.ignore(256, ':');
            inFile >> displayIndividualFitsChoice;
            inFile.ignore(256, ':');
            inFile >> lowerRunHistoIndex;
            inFile.ignore(256, ':');
            inFile >> upperRunHistoIndex;

            fitOptions->SetNumRuns(numRuns);
            fitOptions->SetMultiSourceChoice(true);
            
            element = new ElementFit(BatemanDecaybyActivity, IntegralDecaybyActivity, batemanFitFunctions, integralFitFunctions, paraVals, fitOptions);
            Run* elementRun = new Run(element); 
            
            elementRun->RunTotalBatemanRunsNoChange(0);
            elementRun->RunTotalIntegralRunsNoChange(0);
            if(singleElementDataChoice == 2)
            {
                elementRun->RunSingleBatemanRunsNoChange(0);
                elementRun->RunSingleIntegralRunsNoChange(0);
            }

            //choice for displaying data via histogram
            if(graphOrHistoChoice == 1)
            {
                //calls functions to create and fill histograms with data
                elementRun->CreateTotalBatemanMultiRunHistos();
                elementRun->CreateTotalIntegralMultiRunHistos();
                elementRun->FillTotalBatemanMultiRunHistos();
                elementRun->FillTotalIntegralMultiRunHistos();
                if(singleElementDataChoice == 2)
                {
                    elementRun->CreateSingleBatemanMultiRunHistos();
                    elementRun->CreateSingleIntegralMultiRunHistos();
                    elementRun->FillSingleBatemanMultiRunHistos();
                    elementRun->FillSingleIntegralMultiRunHistos();
                }
                //creates canvases to display histos and displays them
                TCanvas** runResultCanvases = new TCanvas* [numElements];
                for(int i = 0;i < numElements; i++)
                {
                    if(singleElementDataChoice != 2)
                    {
                        runResultCanvases[i] = new TCanvas((elementNames[i] + "ResultHisto").c_str(), (elementNames[i] + "ResultHisto").c_str(), 1100, 550);
                        runResultCanvases[i]->Divide(2,1,.02,.02);
                    }else{
                        runResultCanvases[i] = new TCanvas((elementNames[i] + "ResultHisto").c_str(), (elementNames[i] + "ResultHisto").c_str(), 1100, 1100);
                        runResultCanvases[i]->Divide(2,2,.02,.02);
                    }
                }
                elementRun->DisplayTotalBatemanFitValuesHistos(runResultCanvases);
                elementRun->DisplayTotalIntegralFitValuesHistos(runResultCanvases);
                if(singleElementDataChoice == 2)
                {
                    elementRun->DisplaySingleBatemanFitValuesHistos(runResultCanvases);
                    elementRun->DisplaySingleIntegralFitValuesHistos(runResultCanvases);
                }

                delete [] runResultCanvases;
            //choice of displaying data via graph
            }else{
                elementRun->GenTotalBatemanGraphsNoChange();              
                elementRun->GenTotalIntegralGraphsNoChange();
                if(singleElementDataChoice == 2)
                {
                    elementRun->GenSingleBatemanGraphsNoChange();
                    elementRun->GenSingleIntegralGraphsNoChange();
                }

                TCanvas** runResultCanvases = new TCanvas* [numElements];
                for(int i = 0; i < numElements; i++)
                {
                    if(singleElementDataChoice != 2)
                    {
                        runResultCanvases[i] = new TCanvas((elementNames[i] + "ResultHisto").c_str(), (elementNames[i] + "ResultHisto").c_str(), 1100, 550);
                        runResultCanvases[i]->Divide(2,1,.02,.02);
                    }else{
                        runResultCanvases[i] = new TCanvas((elementNames[i] + "ResultHisto").c_str(), (elementNames[i] + "ResultHisto").c_str(), 1100, 1100);
                        runResultCanvases[i]->Divide(2,2,.02,.02);
                    }
                }
                elementRun->DisplayTotalBatemanFitValuesGraphs(runResultCanvases);
                elementRun->DisplayTotalIntegralFitValuesGraphs(runResultCanvases);
                if(singleElementDataChoice == 2)
                {
                    elementRun->DisplaySingleBatemanFitValuesGraphs(runResultCanvases);
                    elementRun->DisplaySingleIntegralFitValuesGraphs(runResultCanvases);
                }

                delete [] runResultCanvases;
            }

        //case for displaying the individual fits for runs
        if(displayIndividualFitsChoice == 1)
        {
            CycleCanvasHolder* batemanCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, "Total Bateman Fit");
            CycleCanvasHolder* integralCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, "Total Integral Fit");
            
            element->DrawTotalBatemanIndividualHistos(batemanCanvases, lowerRunHistoIndex, upperRunHistoIndex, 0, 0);
            element->DrawTotalIntegralIndividualHistos(integralCanvases, lowerRunHistoIndex, upperRunHistoIndex, 0, 0);

            delete batemanCanvases;
            delete integralCanvases;

            if(singleElementDataChoice == 2)
            {
                SingleCycleCanvasHolder* singleBatemanCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, numElements, elementNames, "Single bateman Fit");
                SingleCycleCanvasHolder* singleIntegralCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, numElements, elementNames, "Single Integral Fit");

                element->DrawSingleBatemanIndividualHistos(singleBatemanCanvases, lowerRunHistoIndex, upperRunHistoIndex, 0, 0);
                element->DrawSingleIntegralIndividualHistos(singleIntegralCanvases, lowerRunHistoIndex, upperRunHistoIndex, 0, 0);
                
                delete singleBatemanCanvases;
                delete singleIntegralCanvases;
            }
        }

        inFile.close();
        delete element;
        delete elementRun;
        break;
        }

        //multiple cycles of histogram
        case 3:
        {
            Int_t numRuns, numCycles, singleSourceHistoChoice, eventChangeChoice, eventChangeType, eventChangeIncrement, timeFitChoice,
                  timeShiftType, displayFitAveragesChoice, lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex;
            Double_t binTimeFitInc, runMeanDifference;

            inFile.open("simulatedMultiCycle.txt");
            inFile.ignore(256,':');
            inFile >> numRuns;
            inFile.ignore(256,':');
            inFile >> numCycles;
            inFile.ignore(256,':');
            inFile >> timeFitChoice;
            inFile.ignore(256,':');
            inFile >> timeShiftType;
            inFile.ignore(256,':');
            inFile >> binTimeFitInc;
            inFile.ignore(256,':');
            inFile >> rebinChoice;
            inFile.ignore(256,':');
            inFile >> rebinBinInc;
            inFile.ignore(256,':');
            inFile >> eventChangeChoice;
            inFile.ignore(256,':');
            inFile >> eventChangeType;
            inFile.ignore(256,':');
            inFile >> eventChangeIncrement;
            inFile.ignore(256,':');
            inFile >> singleSourceHistoChoice;
            inFile.ignore(256,':');
            inFile >> runMeanDifference;
            inFile.ignore(256,':');
            inFile >> displayIndividualFitsChoice;
            inFile.ignore(256,':');
            inFile >> displayFitAveragesChoice;
            inFile.ignore(256,':');
            inFile >> lowerRunHistoIndex;
            inFile.ignore(256,':');
            inFile >> upperRunHistoIndex;
            inFile.ignore(256,':');
            inFile >> lowerCycleHistoIndex;
            inFile.ignore(256,':');
            inFile >> upperCycleHistoIndex;

            fitOptions->SetNumRuns(numRuns);
            fitOptions->SetNumCycles(numCycles);
            fitOptions->SetTimeShiftType(timeShiftType);
            fitOptions->SetTimeFitBinInc(binTimeFitInc);
            fitOptions->SetEventChangeType(eventChangeType);
            fitOptions->SetEventChangeFactor(eventChangeIncrement);
            if(singleSourceHistoChoice == 2)
            {
                fitOptions->SetMultiSourceChoice(true);
            }
            if(runMeanDifference == 1)
            {
                fitOptions->SetRunMeanDifference(true);
            }
            if(rebinChoice == 1)
            {
                fitOptions->SetRebinChoice(true);
            }
            if(eventChangeChoice == 1)
            {
                fitOptions->SetEventNumChangeChoice(true);
            }
            if(timeFitChoice == 1)
            {
                fitOptions->SetTimeFitChoice(true);
            }
            if(displayFitAveragesChoice == 1)
            {
                fitOptions->SetDisplayFitAverages(true);
            }
            fitOptions->SetRebinBinInc(rebinBinInc);

            element = new ElementFit(BatemanDecaybyActivity, IntegralDecaybyActivity, batemanFitFunctions, integralFitFunctions, paraVals, fitOptions);
            Run* elementRunsCycle = new Run(element);
            Cycle* cycle = new Cycle(elementRunsCycle, element);

            TCanvas** canvasArr;

            //time fit change single histo source
            if(fitOptions->GetMultiSourceChoice() == false && fitOptions->GetTimeFitChoice() == true)
            {
                //mean difference
                if(runMeanDifference == 1)
                {   //runs the cycles
                    cycle->RunTotalBatemanCyclesSingleGen();
                    cycle->RunTotalIntegralCyclesSingleGen();
                    //calculates the difference between the bateman and integral halflife fit means for each cycle
                    cycle->GenTotalMeanDifference();
                    //puts the differences in a graph
                    cycle->GenTotalMeanDifferenceGraphsTimeDifference();
                    if(singleElementDataChoice == 2)
                    {   //runs the cycles
                        cycle->RunSingleBatemanCyclesSingleGen();
                        cycle->RunSingleIntegralCyclesSingleGen();
                        //calculates the difference between the bateman and integral halflife fit means for each cycle
                        cycle->GenSingleMeanDifference();
                        //puts the differences in a graph
                        cycle->GenSingleMeanDifferenceGraphsTimeDifference();
                    }

                    canvasArr = new TCanvas* [numElements];
                    for(int i = 0; i < numElements; i++)
                    {
                        if(singleElementDataChoice != 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Single Source Mean Difference").c_str(), (elementNames[i] + " Single Source Mean Difference").c_str(), 500, 500);
                            canvasArr[i]->Divide(1,1,.02,.02);
                        }else if(singleElementDataChoice == 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Single Source Mean Difference").c_str(), (elementNames[i] + " Single Source Mean Difference").c_str(), 1100, 500);
                            canvasArr[i]->Divide(2,1,.02,.02);
                        }
                    }
                    //displays halflife mean difference data
                    cycle->DisplayTotalMeanDifferenceGraphs(canvasArr);
                    if(singleElementDataChoice == 2)
                    {
                        cycle->DisplaySingleMeanDifferenceGraphs(canvasArr);
                    }
                    delete [] canvasArr;
                //seperate mean
                }else{
                    //runs the cycles
                    cycle->RunTotalBatemanCyclesSingleGen();
                    cycle->RunTotalIntegralCyclesSingleGen();
                    //puts the halflife mean fit data in respective graph
                    cycle->GenTotalBatemanMeanGraphsTimeDifference();
                    cycle->GenTotalIntegralMeanGraphsTimeDifference();

                    if(singleElementDataChoice == 2)
                    {   //runs the cycles
                        cycle->RunSingleBatemanCyclesSingleGen();
                        cycle->RunSingleIntegralCyclesSingleGen();
                        //puts the halflife mean fit data in respective graph
                        cycle->GenSingleBatemanMeanGraphsTimeDifference();
                        cycle->GenSingleIntegralMeanGraphsTimeDifference();
                    }

                    canvasArr = new TCanvas* [numElements];
                    for(int i = 0; i < numElements; i++)
                    {
                        if(singleElementDataChoice != 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Single Source Seperate Mean").c_str(), (elementNames[i] + " Single Source Seperate Mean").c_str(), 1100, 500);
                            canvasArr[i]->Divide(2,1,.02,.02);
                        }else if(singleElementDataChoice == 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Single Source Seperate Mean").c_str(), (elementNames[i] + " Single Source Seperate Mean").c_str(), 1100, 1100);
                            canvasArr[i]->Divide(2,2,.02,.02);
                        }
                    }
                    //displays halflife fit data
                    cycle->DisplayTotalBatemanMeanGraphs(canvasArr);
                    cycle->DisplayTotalIntegralMeanGraphs(canvasArr);
                    if(singleElementDataChoice == 2)
                    {
                        cycle->DisplaySingleBatemanMeanGraphs(canvasArr);
                        cycle->DisplaySingleIntegralMeanGraphs(canvasArr);
                    }
                    delete [] canvasArr;
                }
            //time fit change multi histo source
            }else if(fitOptions->GetMultiSourceChoice() == true && fitOptions->GetTimeFitChoice() == true)
            {
                //mean difference
                if(runMeanDifference == 1)
                {   //runs the cycles and calculates the difference between the bateman and integral halflife fit means for each cycle
                    cycle->RunTotalDifferenceMeanCycles();
                    //takes the difference mean halflife data and puts it into a graph
                    cycle->GenTotalMeanDifferenceGraphsTimeDifference();

                    if(singleElementDataChoice == 2)
                    {   //runs the cycles and calculates the difference between the bateman and integral halflife fit means for each cycle
                        cycle->RunSingleDifferenceMeanCycles();
                        //takes the difference mean halflife data and puts it into a graph
                        cycle->GenSingleMeanDifferenceGraphsTimeDifference();
                    }

                    canvasArr = new TCanvas* [numElements];
                    for(int i = 0; i < numElements; i++)
                    {
                        if(singleElementDataChoice != 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Multi Source Mean Difference").c_str(), (elementNames[i] + " Multi Source Mean Difference").c_str(), 500, 500);
                            canvasArr[i]->Divide(1,1,.02,.02);
                        }else if(singleElementDataChoice == 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Multi Source Mean Difference").c_str(), (elementNames[i] + " Multi Source Mean Difference").c_str(), 1100, 500);
                            canvasArr[i]->Divide(2,1,.02,.02);
                        }
                    }
                    cycle->DisplayTotalMeanDifferenceGraphs(canvasArr);
                    if(singleElementDataChoice == 2)
                    {
                        cycle->DisplaySingleMeanDifferenceGraphs(canvasArr);
                    }
                    delete [] canvasArr;
                //seperate mean
                }else if(runMeanDifference == 2)
                {
                    cycle->RunTotalBatemanCycles();
                    cycle->RunTotalIntegralCycles();

                    cycle->GenTotalBatemanMeanGraphsTimeDifference();
                    cycle->GenTotalIntegralMeanGraphsTimeDifference();
                    if(singleElementDataChoice == 2)
                    {
                        cycle->RunSingleBatemanCycles();
                        cycle->RunSingleIntegralCycles();

                        cycle->GenSingleBatemanMeanGraphsTimeDifference();
                        cycle->GenSingleIntegralMeanGraphsTimeDifference();
                    }

                    canvasArr = new TCanvas* [numElements];
                    for(int i = 0; i < numElements; i++)
                    {
                        if(singleElementDataChoice != 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Multi Source Seperate Mean").c_str(), (elementNames[i] + " Multi Source Seperate Mean").c_str(), 1100, 500);
                            canvasArr[i]->Divide(2,1,.02,.02);
                        }else if(singleElementDataChoice == 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Multi Source Seperate Mean").c_str(), (elementNames[i] + " Multi Source Seperate Mean").c_str(), 1100, 1100);
                            canvasArr[i]->Divide(2,2,.02,.02);
                        }
                    }
                    
                    cycle->DisplayTotalBatemanMeanGraphs(canvasArr);
                    cycle->DisplayTotalIntegralMeanGraphs(canvasArr);
                    if(singleElementDataChoice == 2)
                    {
                        cycle->DisplaySingleBatemanMeanGraphs(canvasArr);
                        cycle->DisplaySingleIntegralMeanGraphs(canvasArr);
                    }
                    delete [] canvasArr;
                }
            //rebinning option
            }else if(fitOptions->GetRebinChoice() == true)
            {
                //mean difference
                if(runMeanDifference == 1)
                {   //runs cycles
                    cycle->RunTotalBatemanCycles();
                    cycle->RunTotalIntegralCycles();
                    //takes difference of bateman and integral mean fit half life values
                    cycle->GenTotalMeanDifference();
                    //puts difference mean half life values and puts in a graph
                    cycle->GenTotalMeanDifferenceGraphsRebin();

                    if(singleElementDataChoice == 2)
                    {   //runs cycles
                        cycle->RunSingleBatemanCycles();
                        cycle->RunSingleIntegralCycles();
                        //takes difference of bateman and integral mean fit half life values
                        cycle->GenSingleMeanDifference();
                        //puts difference mean half life values and puts in a graph
                        cycle->GenSingleMeanDifferenceGraphsRebin();
                    }

                    canvasArr = new TCanvas* [numElements];
                    for(int i = 0; i < numElements; i++)
                    {
                        if(singleElementDataChoice != 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Rebin Difference Mean").c_str(), (elementNames[i] + " Rebin Difference Mean").c_str(), 500, 500);
                            canvasArr[i]->Divide(1,1,.02,.02);
                        }else if(singleElementDataChoice == 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Rebin Difference Mean").c_str(), (elementNames[i] + " Rebin Difference Mean").c_str(), 1100, 500);
                            canvasArr[i]->Divide(2,1,.02,.02);
                        }
                    }

                    cycle->DisplayTotalMeanDifferenceGraphs(canvasArr);
                    if(singleElementDataChoice == 2)
                    {
                        cycle->DisplaySingleMeanDifferenceGraphs(canvasArr);
                    }
                    delete [] canvasArr;
                //seperate mean
                }else if(runMeanDifference == 2)
                {   //runs cycles
                    cycle->RunTotalBatemanCycles();
                    cycle->RunTotalIntegralCycles();
                    //puts bateman and integral mean fit data in a graph
                    cycle->GenTotalBatemanMeanGraphsRebin();
                    cycle->GenTotalIntegralMeanGraphsRebin();
                    if(singleElementDataChoice == 2)
                    {   //runs cycles
                        cycle->RunSingleBatemanCycles();
                        cycle->RunSingleIntegralCycles();
                        //puts bateman and integral mean fit data in a graph
                        cycle->GenSingleBatemanMeanGraphsRebin();
                        cycle->GenSingleIntegralMeanGraphsRebin();
                    }

                    canvasArr = new TCanvas* [numElements];
                    for(int i = 0; i < numElements; i++)
                    {
                        if(singleElementDataChoice != 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Rebin Seperate Mean").c_str(), (elementNames[i] + " Rebin Seperate Mean").c_str(), 1100, 500);
                            canvasArr[i]->Divide(2,1,.02,.02);
                        }else if(singleElementDataChoice == 2)
                        {
                            canvasArr[i] = new TCanvas((elementNames[i] + " Rebin Seperate Mean").c_str(), (elementNames[i] + " Rebin Seperate Mean").c_str(), 1100, 1100);
                            canvasArr[i]->Divide(2,2,.02,.02);
                        }
                    }
                    //displays data
                    cycle->DisplayTotalBatemanMeanGraphs(canvasArr);
                    cycle->DisplayTotalIntegralMeanGraphs(canvasArr);
                    if(singleElementDataChoice == 2)
                    {
                        cycle->DisplaySingleBatemanMeanGraphs(canvasArr);
                        cycle->DisplaySingleIntegralMeanGraphs(canvasArr);
                    }
                    delete [] canvasArr;
                }
            //Changing events between cycles
            }else if(fitOptions->GetNumEventChangeChoice() == true)
            {
                Double_t* eventNums = fitOptions->GetEventNumArr();
                for(int i = 0; i < numCycles; i++)
                {
                    cout << eventNums[i] << endl;
                }

            }

            //case for displaying the individual fits for runs
            if(displayIndividualFitsChoice == 1)
            {
                CycleCanvasHolder* batemanCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, "Total Bateman Fit");
                CycleCanvasHolder* integralCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, "Total Integral Fit");
                
                element->DrawTotalBatemanIndividualHistos(batemanCanvases, lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex);
                element->DrawTotalIntegralIndividualHistos(integralCanvases, lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex);

                delete batemanCanvases;
                delete integralCanvases;

                if(singleElementDataChoice == 2)
                {
                    SingleCycleCanvasHolder* singleBatemanCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, numElements, elementNames, "Single bateman Fit");
                    SingleCycleCanvasHolder* singleIntegralCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, numElements, elementNames, "Single Integral Fit");

                    element->DrawSingleBatemanIndividualHistos(singleBatemanCanvases, lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex);
                    element->DrawSingleIntegralIndividualHistos(singleIntegralCanvases, lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex);
                    
                    delete singleBatemanCanvases;
                    delete singleIntegralCanvases;
                }
            }

            delete element;
            delete elementRunsCycle;
            delete cycle;
            inFile.close();
        }
    }
    //delete dyanmically allocated data
    for(int i = 0; i < numElements; i++)
    {
        delete paraVals[i];
    }
    delete [] paraVals;
    delete [] elementNames;
    delete fitOptions;
}

Double_t BatemanDecaybyActivity(Double_t *x, Double_t *par)
{
    Double_t hold = 0.0;
    for(int i = 0; i < numElements; i++)
    {
        hold += batemanFitFunctions[i](x, par);
    }
    return hold;
}

Double_t IntegralDecaybyActivity(Double_t *x, Double_t *par)
{
    Double_t hold = 0.0;
    for(int i = 0; i < numElements; i++)
    {
        hold += integralFitFunctions[i](x, par);
    }
    return hold;
}

//used to create the fit function dynamically
//just have faith the overly complicated conditionals are correct
void CreateFitFunctions(string* elementNames)
{
    ofstream outFile;
    outFile.open("FitFunction.h");
    string* N0Names = new string [numElements];
    string* lambdaNames = new string [numElements];
    string* batemanFunctionNames = new string[numElements];
    string* integralFunctionNames = new string[numElements];

    for(int i = 0; i < numElements; i++)
    {
        N0Names[i] = elementNames[i] + "0";
        lambdaNames[i] = "lambda" + elementNames[i];
        batemanFunctionNames[i] = elementNames[i] + "DecayByActivity";
        integralFunctionNames[i] = elementNames[i] + "DecayByActivityIntegral";
    }

    //class definition
    outFile << "#ifndef FITFUNCTION_H" << endl << "#define FITFUNCTION_H" << endl << endl;
    outFile << "#include \"TMath.h\"" << endl << endl;
    outFile << "using namespace std;" << endl << endl;

    outFile << "class FitFunction{" << endl;
    outFile << "public:" << endl;
    outFile << "\ttypedef Double_t (*decayFunction)(Double_t *x, Double_t *par);" << endl;
    for(int i = 0; i < numElements; i++)
    {
        outFile << "\tstatic Double_t " << batemanFunctionNames[i] << "(Double_t *x, Double_t* par);" << endl;
        outFile << "\tstatic Double_t " << integralFunctionNames[i] << "(Double_t *x, Double_t* par);" << endl;
    }
    outFile << "\tFitFunction(Int_t numElements);" << endl;
    outFile << "\tdecayFunction* GetBatemanFitFunctions(){return batemanFitFunctions;}" << endl;
    outFile << "\tdecayFunction* GetIntegralFitFunctions(){return integralFitFunctions;}" << endl;

    outFile << "private:" << endl;
    outFile << "\tInt_t numElements;" << endl;
    outFile << "\tdecayFunction* batemanFitFunctions;" << endl;
    outFile << "\tdecayFunction* integralFitFunctions;" << endl;
    outFile << "};" << endl << endl;

    //constructor
    outFile << "FitFunction::FitFunction(Int_t numElements)" << endl;
    outFile << "{" << endl;
    outFile << "\tthis->numElements = numElements;" << endl;
    outFile << "\tbatemanFitFunctions = new decayFunction [numElements];" << endl;
    outFile << "\tintegralFitFunctions = new decayFunction [numElements];" << endl;
    for(int i = 0; i < numElements; i++)
    {
        outFile << "\tbatemanFitFunctions[" << i << "] = " << batemanFunctionNames[i] << ";" << endl;
        outFile << "\tintegralFitFunctions[" << i << "] = " << integralFunctionNames[i] << ";" << endl; 
    }
    outFile << "}" << endl << endl;

    //create functions
    for(int elementIndex = 0; elementIndex < numElements; elementIndex++)
    {
        //FOR BATEMAN FUNCTION
        outFile << "Double_t FitFunction::" << batemanFunctionNames[elementIndex] << "(Double_t *x, Double_t *par)" << endl;
        outFile << "{" << endl;

        //variable declaration
        outFile << "\tFloat_t timeVar = x[0];" << endl;
        for(int elementSubIndex = 0; elementSubIndex <= elementIndex; elementSubIndex++)
        {
            outFile << "\tDouble_t " << N0Names[elementSubIndex] << " = par["
            << (elementSubIndex * 2) << "];" << endl;
            outFile << "\tDouble_t " << lambdaNames[elementSubIndex] << " = par["
            << ((elementSubIndex*2) + 1) << "];" << endl;
        }
        outFile << endl;

        //fit function creation
        outFile << "\tDouble_t f = (" << N0Names[elementIndex] << " * " << lambdaNames[elementIndex] << " * (TMath::Exp(-" 
        << lambdaNames[elementIndex] << " * timeVar)));" << endl << endl;

        for(int m = 0; m < elementIndex; m++)
        {
            for(int k = m; k <= elementIndex; k++)
            {
                outFile << "\tf += (" << N0Names[m] << " * " << lambdaNames[elementIndex];
                for(int q = m; q < elementIndex; q++)
                {
                    outFile << " * " << lambdaNames[q];
                }
                outFile << " * ((TMath::Exp(-" << lambdaNames[k] << " * timeVar)) / (";
                if(m != k)
                {
                    outFile << "(" << lambdaNames[m] << "-" << lambdaNames[k] << ")";
                }
                for(int j = m+1; j <= elementIndex; j++)
                {
                    if((j != k && k != m) || (j != k && j != (m+1) && k == m))
                    {
                        outFile << "*(" << lambdaNames[j] << "-" << lambdaNames[k] << ")";
                    }
                    if(j != k && k == m && j == (m+1))
                    {
                        outFile << "(" << lambdaNames[j] << "-" << lambdaNames[k] << ")";
                    }
                }
                outFile << ")));" << endl;
            }
            outFile << endl;
        }
        outFile << "\treturn f;" << endl << "}" << endl << endl;


        //FOR INTEGRAL FUNCTION
        outFile << "Double_t FitFunction::" << integralFunctionNames[elementIndex] << "(Double_t *x, Double_t *par)" << endl;
        outFile << "{" << endl;

        //variable declaration
        outFile << "\tFloat_t timeVar = x[0];" << endl;
        for(int elementSubIndex = 0; elementSubIndex <= elementIndex; elementSubIndex++)
        {
            outFile << "\tDouble_t " << N0Names[elementSubIndex] << " = par["
            << (elementSubIndex * 2) << "];" << endl;
            outFile << "\tDouble_t " << lambdaNames[elementSubIndex] << " = par["
            << ((elementSubIndex*2) + 1) << "];" << endl;
        }
        outFile << endl;

        //function creation
        outFile << "\tDouble_t f = " << N0Names[elementIndex] << " * (1.0 - TMath::Exp(-" << lambdaNames[elementIndex] << " * timeVar));" << endl << endl;

        for(int m = 0; m < elementIndex; m++)
        {
            for(int k = m; k <= elementIndex; k++)
            {
                outFile << "\tf += " << N0Names[m]; 
                for(int q = m; q <= elementIndex; q++)
                {
                    if(q != k)
                    {
                        outFile << " * " << lambdaNames[q];
                    }
                }
                outFile << " * (1.0 - TMath::Exp(-" << lambdaNames[k] << " * timeVar)) / (";
                if(m != k)
                {
                    outFile << "(" << lambdaNames[m] << "-" << lambdaNames[k] << ")";
                }
                for(int j = m+1; j <= elementIndex; j++)
                {
                    if((j != k && k != m) || (j != k && j != (m+1) && k == m))
                    {
                        outFile << "*(" << lambdaNames[j] << "-" << lambdaNames[k] << ")";
                    }
                    if(j != k && k == m && j == (m+1))
                    {
                        outFile << "(" << lambdaNames[j] << "-" << lambdaNames[k] << ")";
                    }
                }
                outFile << ");" << endl;
            }
            outFile << endl;
        }
        outFile << "\treturn f;" << endl << "}" << endl << endl;
    }

    outFile << "#endif";
    delete [] N0Names;
    delete [] lambdaNames;
    delete [] batemanFunctionNames;
    delete [] integralFitFunctions;
}
