#include <fstream>
#include "TMath.h"
#include "TCanvas.h"
#include "ElementFit.h"
#include "TGraph.h"
#include "TH1D.h"
#include "Run.h"
#include "Cycle.h"
#include "TGraphErrors.h"
#include "ParameterValue.h"
#include "TF1.h"
#include "CycleCanvasHolder.h"
#include "SingleCycleCanvasHolder.h"

using namespace std;

Double_t BADecaybyActivityIntegral(Double_t *x, Double_t *par);
Double_t CSDecaybyActivityIntegral(Double_t *x, Double_t *par);
Double_t LADecaybyActivityIntegral(Double_t *x, Double_t *par);
Double_t CSDecaybyActivity(Double_t* x, Double_t *par);
Double_t BADecaybyActivity(Double_t *x, Double_t *par);
Double_t LADecaybyActivity(Double_t *x, Double_t *par);
Double_t IntegralDecaybyActivity(Double_t *x, Double_t *par);
Double_t RegularDecaybyActivity(Double_t *x, Double_t *par);

typedef Double_t (*decayFunction)(Double_t *x, Double_t *par);
decayFunction* regularFitFunctions;
decayFunction* integralFitFunctions;
int numElements;

void fitMultiple();

int main()
{
    fitMultiple();
    return 0;
}

//RUN(N)- process that happens N times in which histograms of the decay is generated and the fit parameters are extracted and put into a histogram
//CYCLE(M)- process that happens M time in which the mean values for the histograms of the RUN process are in some way compared and then graphed
void fitMultiple()
{
    Int_t eventDecrement, numBins, rebinDifference;
    Double_t timeRun, valueHolder, leaveOutStartBinNumber, leaveOutEndBinNumber;
    int histoSimulateChoice, rangeValueChoice, displayIndividualFitsChoice, rebinChoice;
    ElementFit* element;
    ifstream inFile;
    
    //open first option file
    inFile.open("elementInput.txt");
    inFile.ignore(256,':');
    inFile >> numElements;

    //fitFunctions array passed in and used in element object
    Int_t numFitFunction = numElements*2;
    decayFunction* fitFunctions = new decayFunction[numFitFunction];
    fitFunctions[0] = CSDecaybyActivity;
    fitFunctions[1] = CSDecaybyActivityIntegral;
    fitFunctions[2] = BADecaybyActivity;
    fitFunctions[3] = BADecaybyActivityIntegral;
    fitFunctions[4] = LADecaybyActivity;
    fitFunctions[5] = LADecaybyActivityIntegral;
    
    //arrays used to create the full decay function in a modular fassion in main
    regularFitFunctions = new decayFunction [numElements];
    regularFitFunctions[0] = CSDecaybyActivity;
    regularFitFunctions[1] = BADecaybyActivity;
    regularFitFunctions[2] = LADecaybyActivity;
    integralFitFunctions = new decayFunction [numElements];
    integralFitFunctions[0] = CSDecaybyActivityIntegral;
    integralFitFunctions[1] = BADecaybyActivityIntegral;
    integralFitFunctions[2] = LADecaybyActivityIntegral;

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
    inFile.close();
    
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
            Int_t events;
            Int_t numBins;
            int fitGenerationChoice;
            int writeToFileChoice;

            inFile.open("simulate.txt");
            //gets events, bins, and timeRun
            inFile.ignore(256,':');
            inFile >> events;
            inFile.ignore(256,':');
            inFile >> numBins;
            inFile.ignore(256,':');
            inFile >> timeRun;
            //gets the fit generation choice
            inFile.ignore(256,':');
            inFile >> fitGenerationChoice;
            inFile.ignore(256,':');
            inFile >> leaveOutStartBinNumber;
            inFile.ignore(256,':');
            inFile >> leaveOutEndBinNumber;
            inFile.close();

            switch(fitGenerationChoice)
            {
                //single run of histogam
                case 1:
                {
                    element = new ElementFit(events, 1, 1, RegularDecaybyActivity, IntegralDecaybyActivity, fitFunctions, numElements, timeRun, numBins, elementNames, paraVals, 2, 2, 0, leaveOutStartBinNumber, leaveOutEndBinNumber);
                    inFile.open("simulated_single_run.txt");
                    element->fitHistos(0, 0);
                    element->displayParameters();
                    //gets write file choice
                    inFile.ignore(256,':');
                    inFile >> writeToFileChoice;
                    //case for writing to file
                    if(writeToFileChoice == 1)
                    {
                        //gets file name
                        string fileName;
                        inFile.ignore(256,';');
                        inFile >> fileName;
                    //case for displaying histogram
                    }else{
                        TCanvas** singleCanvas = new TCanvas* [numElements*2];
                        for(int i = 0; i < numElements; i++)
                        {
                            singleCanvas[(i*2)] = new TCanvas((elementNames[i] + " Single Regular Histo").c_str(), (elementNames[i] + " Single Regular Histo").c_str(), 500, 500);
                            singleCanvas[(i*2)+1] = new TCanvas((elementNames[i] + " Single Integral Histo").c_str(), (elementNames[i] + " Single Integral Histo").c_str(), 500, 500);
                        }
                        TCanvas* regularTotalCanvas = new TCanvas("Total Regular Histo", "Total Regular Histo", 500, 500);
                        TCanvas* integralTotalCanvas = new TCanvas("Total Integral Histo", "Total Integral Histo", 500, 500);
                        element->displaySingleHistos(singleCanvas);
                        element->displayRegularHisto(regularTotalCanvas);
                        element->displayIntegralHisto(integralTotalCanvas);
                        delete [] singleCanvas;
                    }
                    inFile.close();
                    delete element;
                break;
                }
                //multi run of histogram
                case 2:
                {
                    Int_t runs;
                    int graphOrHistoChoice, lowerRunHistoIndex, upperRunHistoIndex;

                    inFile.open("simulated_single_cylce.txt");
                    inFile.ignore(256,':');
                    inFile >> runs;
                    inFile.ignore(256,':');
                    inFile >> eventDecrement;
                    inFile.ignore(256,':');
                    inFile >> writeToFileChoice;
                    inFile.ignore(256,':');
                    inFile >> graphOrHistoChoice;
                    inFile.ignore(256, ':');
                    inFile >> displayIndividualFitsChoice;
                    inFile.ignore(256, ':');
                    inFile >> lowerRunHistoIndex;
                    inFile.ignore(256, ':');
                    inFile >> upperRunHistoIndex;
                    
                    element = new ElementFit(events, runs, 1, RegularDecaybyActivity, IntegralDecaybyActivity, fitFunctions, numElements, timeRun, numBins, elementNames, paraVals, 2, 2, 0, leaveOutStartBinNumber, leaveOutEndBinNumber);
                    Run* elementRun = new Run(runs, eventDecrement, element, elementNames); 
                    
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
                        TGraph** runResultGraphs = elementRun->genGraphsNoChange();
                        TGraph** runResultSingleElementGraphs = elementRun->genGraphsNoChangeSingleElement();

                        TGraph** totalRunResults = new TGraph* [numElements*4];
                        TGraph** totalRunResultErrors = new TGraph* [numElements*4];
                
                
                        //extracts the fit value graphs from the array of graphs runResultGraphs and runResultSingleElementGraphs
                        //(I wanted to put the graphs of single element fit value and total element fit value into one array, but
                        //the graph array runResultGraphs and runResultSingleElementGraphs both had 18 graphs in them with a lot of extra data
                        //so I had to extract the graphs in this really weird way)
                        for(int i = 0; i < numElements; i++)
                        { 
                            totalRunResults[(i*4)] = runResultGraphs[(i*6)+2];
                            totalRunResults[(i*4)+1] = runResultSingleElementGraphs[(i*6)+2];
                            totalRunResults[(i*4)+2] = runResultGraphs[(i*6)+3];
                            totalRunResults[(i*4)+3] = runResultSingleElementGraphs[(i*6)+3];
                        }
                        for(int i = 0; i < numElements; i++)
                        { 
                            totalRunResultErrors[(i*4)] = runResultGraphs[(i*6)];
                            totalRunResultErrors[(i*4)+1] = runResultSingleElementGraphs[(i*6)];
                            totalRunResultErrors[(i*4)+2] = runResultGraphs[(i*6)+1];
                            totalRunResultErrors[(i*4)+3] = runResultSingleElementGraphs[(i*6)+1];
                        }
                        delete [] runResultGraphs;
                        delete [] runResultSingleElementGraphs;

                        //case for writing to file
                        if(writeToFileChoice == 1)
                        {
                            
                            string fileName;
                            inFile.ignore(256,';');
                            inFile >> fileName;
                            
                        //case for displaying results
                        }else{
                            TCanvas** runResultCanvases = new TCanvas* [numElements];
                            TCanvas** runResultErrorCanvases = new TCanvas* [numElements];
                            for(int i = 0; i < numElements; i++)
                            {
                                runResultCanvases[i] = new TCanvas((elementNames[i] + "ResultGraph").c_str(), (elementNames[i] + "ResultGraph").c_str(), 1100, 1100);
                                runResultCanvases[i]->Divide(2,2,.02,.02);
                            }
                            for(int i = 0; i < numElements; i++)
                            {
                                runResultErrorCanvases[i] = new TCanvas((elementNames[i] + "ResultGraphErrors").c_str(), (elementNames[i] + "ResultGraphErrors").c_str(), 1100, 1100);
                                runResultErrorCanvases[i]->Divide(2,2,.02,.02);
                            }
                            elementRun->displayMultiRunResultGraphs(runResultCanvases,totalRunResults);
                            elementRun->displayMultiRunResultGraphs(runResultErrorCanvases, totalRunResultErrors);

                            delete [] runResultErrorCanvases;
                            delete [] runResultCanvases;
                        }
                        delete [] totalRunResults;
                        delete [] totalRunResultErrors;
                    }

                //case for displaying the individual fits for runs
                if(displayIndividualFitsChoice == 1)
                {
                    CycleCanvasHolder* regularCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, "Total Regular Fit");
                    CycleCanvasHolder* integralCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, "Total Integral Fit");
                    SingleCycleCanvasHolder* singleRegularCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, numElements, elementNames, "Single Regular Fit");
                    SingleCycleCanvasHolder* singleIntegralCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, 0, 0, numElements, elementNames, "Single Integral Fit");

                    element->DrawIndividualHistos(regularCanvases, integralCanvases, singleRegularCanvases, singleIntegralCanvases, lowerRunHistoIndex, upperRunHistoIndex, 0, 0);
                    
                    delete regularCanvases;
                    delete integralCanvases;
                    delete singleRegularCanvases;
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
                    Int_t numRuns, numCycles, writeToFileChoice, singleSourceHistoChoice, cycleChangeChoice,
                          lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex;
                    Double_t x_start, x_stop, x_inc, cycleMeanChoice;

                    inFile.open("simulated_multi_cycle.txt");
                    inFile.ignore(256,':');
                    inFile >> numRuns;
                    inFile.ignore(256,':');
                    inFile >> numCycles;
                    inFile.ignore(256,':');
                    inFile >> eventDecrement;
                    inFile.ignore(256,':');
                    inFile >> writeToFileChoice;
                    inFile.ignore(256,':');
                    inFile >> cycleChangeChoice;
                    inFile.ignore(256,':');
                    inFile >> x_start;
                    inFile.ignore(256,':');
                    inFile >> x_stop;
                    inFile.ignore(256,':');
                    inFile >> x_inc;
                    inFile.ignore(256,':');
                    inFile >> singleSourceHistoChoice;
                    inFile.ignore(256,':');
                    inFile >> cycleMeanChoice;
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

                    element = new ElementFit(events, numRuns, numCycles, RegularDecaybyActivity, IntegralDecaybyActivity, fitFunctions, numElements, timeRun, numBins, elementNames, paraVals, singleSourceHistoChoice, rebinChoice, rebinDifference, leaveOutStartBinNumber, leaveOutEndBinNumber);
                    Run* elementRunsCycle = new Run(numRuns, eventDecrement, element, elementNames);
                    Cycle* cycle = new Cycle(numCycles, elementRunsCycle, element, x_start, x_stop, x_inc, cycleChangeChoice);

                    TCanvas** canvasArr;

                    //single histo source
                    if(singleSourceHistoChoice == 1 && rebinChoice == 2)
                    {
                        //mean difference
                        if(cycleMeanChoice == 1)
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
                        if(cycleMeanChoice == 1)
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
                        }else if(cycleMeanChoice == 2)
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
                        if(cycleMeanChoice == 1)
                        {
                            cycle->runDifferenceMeanRebin();
                            cycle->genSeperateMeanGraphsRebin();
                            canvasArr = new TCanvas* [numElements];
                            for(int i = 0; i < numElements; i++)
                            {
                                canvasArr[i] = new TCanvas((elementNames[i] + " Rebin Difference Mean").c_str(), (elementNames[i] + " Rebin Difference Mean").c_str(), 1100, 500);
                                canvasArr[i]->Divide(2,1,.02,.02);
                            }
                            cycle->displayMeanDifferenceGraphs(canvasArr);
                            delete [] canvasArr;
                        //seperate mean
                        }else if(cycleMeanChoice == 2)
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
                        CycleCanvasHolder* regularCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, "Total Regular Fit");
                        CycleCanvasHolder* integralCanvases = new CycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, "Total Integral Fit");
                        SingleCycleCanvasHolder* singleRegularCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, numElements, elementNames, "Single Regular Fit");
                        SingleCycleCanvasHolder* singleIntegralCanvases = new SingleCycleCanvasHolder(lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex, numElements, elementNames, "Single Integral Fit");

                        element->DrawIndividualHistos(regularCanvases, integralCanvases, singleRegularCanvases, singleIntegralCanvases, lowerRunHistoIndex, upperRunHistoIndex, lowerCycleHistoIndex, upperCycleHistoIndex);
                        
                        delete regularCanvases;
                        delete integralCanvases;
                        delete singleRegularCanvases;
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

//formula for CS activity
Double_t CSDecaybyActivityIntegral(Double_t *x, Double_t *par)
{
    Float_t timeVar = x[0];
    Double_t CS0 = par[0];
    Double_t lambdaCS = par[1];

    Double_t f = CS0 * (1 - TMath::Exp(- lambdaCS * timeVar));
    return f;
}

//formula for BA activity
Double_t BADecaybyActivityIntegral(Double_t *x, Double_t *par)
{
    Float_t timeVar = x[0];
    Double_t CS0 = par[0];
    Double_t lambdaCS = par[1];
    Double_t BA0 = par[2];
    Double_t lambdaBA = par[3];

    //inital BA
    Double_t f = BA0*(1- TMath::Exp(-lambdaBA * timeVar));
    //first summation
    f += (((CS0*lambdaBA)/(lambdaBA-lambdaCS))*(1- TMath::Exp(-lambdaCS * timeVar)));
    //second summation
    f += (((CS0*lambdaCS)/(lambdaCS-lambdaBA))*(1- TMath::Exp(-lambdaBA * timeVar)));

    return f;
}

//formla for LA activity
Double_t LADecaybyActivityIntegral(Double_t *x, Double_t *par)
{
    Float_t timeVar = x[0];
    Double_t CS0 = par[0];
    Double_t lambdaCS = par[1]; 
    Double_t BA0 = par[2];
    Double_t lambdaBA = par[3];
    Double_t LA0 = par[4];
    Double_t lambdaLA = par[5];

    //intial LA
    Double_t f = LA0 * ( 1.0 - TMath::Exp( - lambdaLA * timeVar ) );
    
    //intial BA
    f += BA0 * lambdaLA / (lambdaLA - lambdaBA ) * ( 1.0 - TMath::Exp( - lambdaBA * timeVar ) );
    f += BA0 * lambdaBA / (lambdaBA - lambdaLA ) * ( 1.0 - TMath::Exp( - lambdaLA * timeVar ) );

    //first, second, and third part of summation
    f += CS0 * lambdaCS * lambdaBA / (lambdaCS - lambdaLA ) / ( lambdaBA - lambdaLA ) * ( 1.0 - TMath::Exp( - lambdaLA * timeVar ) );
    f += CS0 * lambdaCS * lambdaLA / (lambdaCS - lambdaBA ) / ( lambdaLA - lambdaBA ) * ( 1.0 - TMath::Exp( - lambdaBA * timeVar ) );
    f += CS0 * lambdaBA * lambdaLA / (lambdaBA - lambdaCS ) / ( lambdaLA - lambdaCS ) * ( 1.0 - TMath::Exp( - lambdaCS * timeVar ) );
 
    return f;
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

Double_t CSDecaybyActivity(Double_t* x, Double_t *par)
{
    Float_t timeVar = x[0];
    Double_t CS0 = par[0];
    Double_t lambdaCS = par[1];

    Double_t f = CS0*lambdaCS*TMath::Exp(- lambdaCS* timeVar);
    return f;
}

Double_t BADecaybyActivity(Double_t *x, Double_t *par)
{
    Float_t timeVar = x[0];
    Double_t CS0 = par[0];
    Double_t lambdaCS = par[1];
    Double_t BA0 = par[2];
    Double_t lambdaBA = par[3];

    Double_t f = BA0 * lambdaBA * TMath::Exp(- lambdaBA*timeVar);

    f += CS0 * lambdaCS * lambdaBA * ((TMath::Exp(- lambdaCS*timeVar))/(lambdaBA-lambdaCS));
    f += CS0 * lambdaCS * lambdaBA * ((TMath::Exp(- lambdaBA*timeVar))/(lambdaCS-lambdaBA));
    return f;
}

Double_t LADecaybyActivity(Double_t *x, Double_t *par)
{
    Float_t timeVar = x[0];
    Double_t CS0 = par[0];
    Double_t lambdaCS = par[1];
    Double_t BA0 = par[2];
    Double_t lambdaBA = par[3];
    Double_t LA0 = par[4];
    Double_t lambdaLA = par[5];

    Double_t f = (LA0 * lambdaLA * (TMath::Exp(-lambdaLA * timeVar)));

    f += (BA0 * lambdaBA *lambdaLA*((TMath::Exp(- lambdaBA * timeVar))/(lambdaLA - lambdaBA)));
    f += (BA0 * lambdaBA *lambdaLA*((TMath::Exp(- lambdaLA * timeVar))/(lambdaBA - lambdaLA)));

    f += (CS0 * lambdaCS * lambdaBA * lambdaLA*((TMath::Exp(- lambdaCS * timeVar))/((lambdaBA-lambdaCS)*(lambdaLA-lambdaCS))));
    f += (CS0 * lambdaCS * lambdaBA * lambdaLA*((TMath::Exp(- lambdaBA * timeVar))/((lambdaCS-lambdaBA)*(lambdaLA-lambdaBA))));
    f += (CS0 * lambdaCS * lambdaBA * lambdaLA*((TMath::Exp(- lambdaLA * timeVar))/((lambdaCS-lambdaLA)*(lambdaBA-lambdaLA))));

    return f;
}

Double_t RegularDecaybyActivity(Double_t *x, Double_t *par)
{
    Double_t hold = 0.0;
    for(int i = 0; i < numElements; i++)
    {
        hold += regularFitFunctions[i](x, par);
    }
    return hold;
}
