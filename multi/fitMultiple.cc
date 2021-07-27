#include <fstream>

#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"

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
    Int_t eventDecrement, numBins, rebinDifference, createFitFunctionsChoice;
    Double_t timeRunEnd, valueHolder, leaveOutStartBinNumber, leaveOutEndBinNumber;
    int histoSimulateChoice, rangeValueChoice, displayIndividualFitsChoice, rebinChoice;
    FitOption* fitOptions = new FitOption();
    ElementFit* element;
    ifstream inFile;

    //open first option file
    inFile.open("elementInput.txt");
    inFile.ignore(256,':');
    inFile >> numElements;
    fitOptions->SetNumElements(numElements);

    //fitFunctions array passed in and used in element object
    Int_t numFitFunction = numElements*2;
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
        //get range value choice for init value
        inFile.ignore(256,':');
        inFile >> rangeValueChoice;
        
        if(rangeValueChoice == 1)
        {
            paraVals[i]->setIsValueInitValue(true);
        }else{
            paraVals[i]->setIsValueInitValue(false);
        }
        //gets value for init value
        if(paraVals[i]->getIsValueInitValue())
        {
            inFile.ignore(256,':');
            inFile >> valueHolder;
            paraVals[i]->setInitValue(valueHolder);
        //gets range for init value
        }else{
            inFile.ignore(256,';');
            inFile >> valueHolder;
            paraVals[i]->setLowerRangeInitValue(valueHolder);
            inFile.ignore(256,';');
            inFile >> valueHolder;
            paraVals[i]->setUpperRangeInitValue(valueHolder);
        }
        //get range value choice for half life
        inFile.ignore(256,':');
        inFile >> rangeValueChoice;
        if(rangeValueChoice == 1)
        {
            paraVals[i]->setIsValueHalfLife(true);
            paraVals[i]->setIsValueDecayConst(true);
        }else{
            paraVals[i]->setIsValueHalfLife(false);
            paraVals[i]->setIsValueDecayConst(false);
        }
        //gets half life value
        if(paraVals[i]->getIsValueHalfLife())
        {
            inFile.ignore(256,':');
            inFile >> valueHolder;
            paraVals[i]->setValueHalfLife(valueHolder);
        //gets half life range
        }else{
            inFile.ignore(256,';');
            inFile >> valueHolder;
            paraVals[i]->setLowerRangeHalfLife(valueHolder);
            inFile.ignore(256,';');
            inFile >> valueHolder;
            paraVals[i]->setUpperRangeHalfLife(valueHolder);
        }
    }
    //gets histogram choice
    inFile.ignore(256,':');
    inFile >> histoSimulateChoice;
    inFile.ignore(256,':');
    inFile >> createFitFunctionsChoice;
    inFile.close();

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
    for(int i = 0; i < numElements; i++)
    {
        valueHolder = (log(2)/paraVals[i]->getLowerRangeHalfLife());
        paraVals[i]->setLowerRangeDecayConst(valueHolder);
        valueHolder = (log(2)/paraVals[i]->getUpperRangeHalfLife());
        paraVals[i]->setUpperRangeDecayConst(valueHolder);
        valueHolder = (log(2)/paraVals[i]->getValueHalfLife());
        paraVals[i]->setValueDecayConst(valueHolder);
    }

    switch(histoSimulateChoice)
    {
        //runs case for simulating the histogram
        case 1:
        {
            Int_t events, numBins;
            Double_t binWidth;
            int programExecutionType;

            inFile.open("simulate.txt");
            //gets events, bins, and timeRun
            inFile.ignore(256,':');
            inFile >> events;
            inFile.ignore(256,':');
            inFile >> numBins;
            inFile.ignore(256,':');
            inFile >> binWidth;
            inFile.ignore(256,':');
            inFile >> timeRunEnd;
            //gets the fit generation choice
            inFile.ignore(256,':');
            inFile >> programExecutionType;
            inFile.ignore(256,':');
            inFile >> leaveOutStartBinNumber;
            inFile.ignore(256,':');
            inFile >> leaveOutEndBinNumber;
            inFile.close();

            fitOptions->SetNumEvents(events);
            fitOptions->SetNumBins(numBins);
            fitOptions->SetBinWidth(binWidth);
            fitOptions->SetTimeRunEnd(timeRunEnd);
            fitOptions->SetProgramExecutionType(programExecutionType);
            fitOptions->SetLeaveOutStartBinNumber(leaveOutStartBinNumber);
            fitOptions->SetLeaveOutEndBinNumber(leaveOutEndBinNumber);
            fitOptions->SetElementNames(elementNames);

            switch(programExecutionType)
            {
                //single run of histogam
                case 1:
                {
                    fitOptions->SetMultiSourceChoice(true);
                    element = new ElementFit(BatemanDecaybyActivity, IntegralDecaybyActivity, batemanFitFunctions, integralFitFunctions, paraVals, fitOptions);
                    element->fitHistos(0, 0);
                    element->displayParameters();
                
                    //case for displaying histogram
                    TCanvas** singleCanvas = new TCanvas* [numElements*2];
                    for(int i = 0; i < numElements; i++)
                    {
                        singleCanvas[(i*2)] = new TCanvas((elementNames[i] + " Single Bateman Histo").c_str(), (elementNames[i] + " Single Bateman Histo").c_str(), 500, 500);
                        singleCanvas[(i*2)+1] = new TCanvas((elementNames[i] + " Single Integral Histo").c_str(), (elementNames[i] + " Single Integral Histo").c_str(), 500, 500);
                    }
                    TCanvas* batemanTotalCanvas = new TCanvas("Total Bateman Histo", "Total Bateman Histo", 500, 500);
                    TCanvas* integralTotalCanvas = new TCanvas("Total Integral Histo", "Total Integral Histo", 500, 500);
                    element->displaySingleHistos(singleCanvas);
                    element->displayBatemanHisto(batemanTotalCanvas);
                    element->displayIntegralGraph(integralTotalCanvas);
                    delete [] singleCanvas;
                    delete element;
                break;
                }
                //multi run of histogram
                case 2:
                {
                    Int_t runs;
                    int graphOrHistoChoice, lowerRunHistoIndex, upperRunHistoIndex;

                    inFile.open("simulated_single_cycle.txt");
                    inFile.ignore(256,':');
                    inFile >> runs;
                    inFile.ignore(256,':');
                    inFile >> eventDecrement;
                    inFile.ignore(256,':');
                    inFile >> graphOrHistoChoice;
                    inFile.ignore(256, ':');
                    inFile >> displayIndividualFitsChoice;
                    inFile.ignore(256, ':');
                    inFile >> lowerRunHistoIndex;
                    inFile.ignore(256, ':');
                    inFile >> upperRunHistoIndex;

                    fitOptions->SetNumRuns(runs);
                    fitOptions->SetEventDecrement(eventDecrement);
                    fitOptions->SetMultiSourceChoice(true);
                    
                    element = new ElementFit(BatemanDecaybyActivity, IntegralDecaybyActivity, batemanFitFunctions, integralFitFunctions, paraVals, fitOptions);
                    Run* elementRun = new Run(element); 
                    
                    elementRun->runNoChange(0);

                    //choice for displaying data via histogram
                    if(graphOrHistoChoice == 1)
                    {
                        //calls functions to create and fill histograms with data
                        TH1D** runResultHisto = elementRun->createRunResultHistos();
                        TH1D** runResultHistoSingleElements = elementRun->createRunResultHistosSingleElements();
                        runResultHisto = elementRun->fillRunResultHistos(runResultHisto);
                        runResultHistoSingleElements = elementRun->fillRunResultHistosSingleElement(runResultHistoSingleElements);

                        //put the histograms in these histo arrays to be able to display them right
                        TH1D** totalResultHisto = new TH1D* [numElements*4];
                        for(int i = 0; i < numElements; i++)
                        {
                            totalResultHisto[(i*4)] = runResultHisto[(i*2)];
                            totalResultHisto[(i*4)+1] = runResultHistoSingleElements[(i*2)];
                            totalResultHisto[(i*4)+2] = runResultHisto[(i*2)+1];
                            totalResultHisto[(i*4)+3] = runResultHistoSingleElements[(i*2)+1];
                        }
                        //delete only runResultHisto and runResultHisto arrays because you do not want to delete the individual histograms
                        delete [] runResultHisto;
                        delete [] runResultHistoSingleElements;

                        //creates canvases to display histos and displays them
                        TCanvas** totalRunResultCanvases = new TCanvas* [numElements];
                        for(int i = 0;i < numElements; i++)
                        {
                            totalRunResultCanvases[i] = new TCanvas((elementNames[i] + "ResultHisto").c_str(), (elementNames[i] + "ResultHisto").c_str(), 1100, 1100);
                            totalRunResultCanvases[i]->Divide(2,2,.02,.02);
                        }

                        elementRun->displayMultiRunResultHistos(totalRunResultCanvases, totalResultHisto);

                        delete [] totalRunResultCanvases;
                        delete [] totalResultHisto;

                    //choice of displaying data via graph
                    }else{
                        elementRun->genGraphsNoChange();
                        elementRun->genGraphsNoChangeSingleElement();                
                        
                        TCanvas** runResultCanvases = new TCanvas* [numElements];
                        for(int i = 0; i < numElements; i++)
                        {
                            runResultCanvases[i] = new TCanvas((elementNames[i] + "ResultGraph").c_str(), (elementNames[i] + "ResultGraph").c_str(), 1100, 1100);
                            runResultCanvases[i]->Divide(2,2,.02,.02);
                        }
                        elementRun->displayMultiRunResultGraphs(runResultCanvases);
                        delete [] runResultCanvases;                        
                    }

                //case for displaying the individual fits for runs
                if(displayIndividualFitsChoice == 1)
                {
                    CycleCanvasHolder* batemanCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, "Total Bateman Fit");
                    CycleCanvasHolder* integralCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, "Total Integral Fit");
                    SingleCycleCanvasHolder* singleBatemanCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, numElements, elementNames, "Single bateman Fit");
                    SingleCycleCanvasHolder* singleIntegralCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, numElements, elementNames, "Single Integral Fit");

                    element->DrawIndividualHistos(batemanCanvases, integralCanvases, singleBatemanCanvases, singleIntegralCanvases, lowerRunHistoIndex, upperRunHistoIndex, 0, 0);
                    
                    delete batemanCanvases;
                    delete integralCanvases;
                    delete singleBatemanCanvases;
                    delete singleIntegralCanvases;
                }

                inFile.close();
                delete element;
                delete elementRun;
                break;
                }

                //multiple runs of histogram
                case 3:
                {
                    Int_t numRuns, numCycles, singleSourceHistoChoice, timeShiftType,
                          lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex;
                    Double_t binTimeFitInc, runMeanDifference;

                    inFile.open("simulated_multi_cycle.txt");
                    inFile.ignore(256,':');
                    inFile >> numRuns;
                    inFile.ignore(256,':');
                    inFile >> numCycles;
                    inFile.ignore(256,':');
                    inFile >> eventDecrement;
                    inFile.ignore(256,':');
                    inFile >> timeShiftType;
                    inFile.ignore(256,':');
                    inFile >> binTimeFitInc;
                    inFile.ignore(256,':');
                    inFile >> singleSourceHistoChoice;
                    inFile.ignore(256,':');
                    inFile >> runMeanDifference;
                    inFile.ignore(256,':');
                    inFile >> displayIndividualFitsChoice;
                    inFile.ignore(256,':');
                    inFile >> lowerRunHistoIndex;
                    inFile.ignore(256,':');
                    inFile >> upperRunHistoIndex;
                    inFile.ignore(256,':');
                    inFile >> lowerCycleHistoIndex;
                    inFile.ignore(256,':');
                    inFile >> upperCycleHistoIndex;
                    inFile.ignore(256,':');
                    inFile >> rebinChoice;
                    inFile.ignore(256,':');
                    inFile >> rebinDifference;

                    fitOptions->SetNumRuns(numRuns);
                    fitOptions->SetNumCycles(numCycles);
                    fitOptions->SetEventDecrement(eventDecrement);
                    fitOptions->SetTimeShiftType(timeShiftType);
                    fitOptions->SetTimeFitBinInc(binTimeFitInc);
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
                    fitOptions->SetRebinDifference(rebinDifference);

                    element = new ElementFit(BatemanDecaybyActivity, IntegralDecaybyActivity, batemanFitFunctions, integralFitFunctions, paraVals, fitOptions);
                    Run* elementRunsCycle = new Run(element);
                    Cycle* cycle = new Cycle(elementRunsCycle, element);

                    TCanvas** canvasArr;

                    //single histo source
                    if(singleSourceHistoChoice == 1 && rebinChoice == 2)
                    {
                        //mean difference
                        if(runMeanDifference == 1)
                        {
                            cycle->runSeperateSingleGen();
                            cycle->genSingleMeanDifference();
                            cycle->genMeanDifferenceGraphsTimeDifference();
                            canvasArr = new TCanvas* [numElements];
                            for(int i = 0; i < numElements; i++)
                            {
                                canvasArr[i] = new TCanvas((elementNames[i] + " Single Source Mean Difference").c_str(), (elementNames[i] + " Single Source Mean Difference").c_str(), 1100, 500);
                                canvasArr[i]->Divide(2,1,.02,.02);
                            }
                            cycle->displayMeanDifferenceGraphs(canvasArr);
                            delete [] canvasArr;
                        //seperate mean
                        }else{
                            cycle->runSeperateSingleGen();
                            cycle->genSeperateMeanGraphsTimeDifference();
                            canvasArr = new TCanvas* [numElements];
                            for(int i = 0; i < numElements; i++)
                            {
                                canvasArr[i] = new TCanvas((elementNames[i] + " Single Source Seperate Mean").c_str(), (elementNames[i] + " Single Source Seperate Mean").c_str(), 1100, 1100);
                                canvasArr[i]->Divide(2,2,.02,.02);
                            }
                            cycle->displayMeanSeperateGraphs(canvasArr);
                            delete [] canvasArr;
                        }
                    //multi histo source
                    }else if(singleSourceHistoChoice == 2 && rebinChoice == 2)
                    {
                        //mean difference
                        if(runMeanDifference == 1)
                        {
                            cycle->runDifferenceMeanTimeDifference();
                            cycle->genMeanDifferenceGraphsTimeDifference();
                            canvasArr = new TCanvas* [numElements];
                            for(int i = 0; i < numElements; i++)
                            {
                                canvasArr[i] = new TCanvas((elementNames[i] + " Multi Source Mean Difference").c_str(), (elementNames[i] + " Multi Source Mean Difference").c_str(), 1100, 500);
                                canvasArr[i]->Divide(2,1,.02,.02);
                            }
                            cycle->displayMeanDifferenceGraphs(canvasArr);
                            delete [] canvasArr;
                        //seperate mean
                        }else if(runMeanDifference == 2)
                        {
                            cycle->runSeperateMeanTimeDifference();
                            cycle->genSeperateMeanGraphsTimeDifference();
                            canvasArr = new TCanvas* [numElements];
                            for(int i = 0; i < numElements; i++)
                            {
                                canvasArr[i] = new TCanvas((elementNames[i] + " Multi Source Seperate Mean").c_str(), (elementNames[i] + " Multi Source Seperate Mean").c_str(), 1100, 1100);
                                canvasArr[i]->Divide(2,2,.02,.02);
                            }
                            
                            cycle->displayMeanSeperateGraphs(canvasArr);
                            delete [] canvasArr;
                        }
                    }else if(rebinChoice == 1)
                    {
                        //mean difference
                        if(runMeanDifference == 1)
                        {
                            cycle->runDifferenceMeanRebin();
                            cycle->genMeanDifferenceGraphsRebin();
                            canvasArr = new TCanvas* [numElements];
                            for(int i = 0; i < numElements; i++)
                            {
                                canvasArr[i] = new TCanvas((elementNames[i] + " Rebin Difference Mean").c_str(), (elementNames[i] + " Rebin Difference Mean").c_str(), 1100, 500);
                                canvasArr[i]->Divide(2,1,.02,.02);
                            }
                            cycle->displayMeanDifferenceGraphs(canvasArr);
                            delete [] canvasArr;
                        //seperate mean
                        }else if(runMeanDifference == 2)
                        {
                            cycle->runSeperateMeanRebin();
                            cycle->genSeperateMeanGraphsRebin();
                            canvasArr = new TCanvas* [numElements];
                            for(int i = 0; i < numElements; i++)
                            {
                                canvasArr[i] = new TCanvas((elementNames[i] + " Rebin Seperate Mean").c_str(), (elementNames[i] + " Rebin Seperate Mean").c_str(), 1100, 500);
                                canvasArr[i]->Divide(2,2,.02,.02);
                            }
                            cycle->displayMeanSeperateGraphs(canvasArr);

                            delete [] canvasArr;
                        }
                    }

                    //case for displaying the individual fits for runs
                    if(displayIndividualFitsChoice == 1)
                    {
                        CycleCanvasHolder* batemanCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, "Total bateman Fit");
                        CycleCanvasHolder* integralCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, "Total Integral Fit");
                        SingleCycleCanvasHolder* singleBatemanCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, numElements, elementNames, "Single bateman Fit");
                        SingleCycleCanvasHolder* singleIntegralCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, numElements, elementNames, "Single Integral Fit");

                        element->DrawIndividualHistos(batemanCanvases, integralCanvases, singleBatemanCanvases, singleIntegralCanvases, lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex);
                        
                        delete batemanCanvases;
                        delete integralCanvases;
                        delete singleBatemanCanvases;
                        delete singleIntegralCanvases;
                    }

                    delete element;
                    delete elementRunsCycle;
                    delete cycle;
                    inFile.close();
                }
            }
            break;
        }
    }
    //delete dyanmically allocated data
    for(int i = 0; i < numElements; i++)
    {
        delete paraVals[i];
    }
    delete [] paraVals;
    delete [] elementNames;
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
