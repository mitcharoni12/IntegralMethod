/*
    PURPOSE: used to store the graphs for each cycle of the integral histograms
*/
#ifndef CYCLEGRAPHHOLDER_H
#define CYCLEGRAPHHOLDER_H

using namespace std;

#include "TGraph.h"
#include "RunGraphHolder.h"

class CycleGraphHolder{
private:
    Int_t numCycles;
    RunGraphHolder** graphArr;
public:
    CycleGraphHolder(Int_t numCycles, Int_t numRuns, string graphName, Int_t* numPoints);
    ~CycleGraphHolder();
    TGraph* GetAGraph(Int_t cycleIndex, Int_t runIndex);
    void SetAGraph(Int_t cycleIndex, Int_t runIndex, TGraph* graph);
};

CycleGraphHolder::CycleGraphHolder(Int_t numCycles, Int_t numRuns, string graphName, Int_t* numPoints)
{
    this->numCycles = numCycles;

    graphArr = new RunGraphHolder* [numCycles];
    for(int i = 0; i < numCycles; i++)
    {
        graphArr[i] = new RunGraphHolder(i, numRuns, graphName, numPoints[i]);
    }
}

CycleGraphHolder::~CycleGraphHolder()
{
    for(int i = 0; i < numCycles; i++)
    {
        delete graphArr[i];
    }
    delete [] graphArr;
}

TGraph* CycleGraphHolder::GetAGraph(Int_t cycleIndex, Int_t runIndex)
{
    return graphArr[cycleIndex]->GetAGraph(runIndex);
}

void CycleGraphHolder::SetAGraph(Int_t cycleIndex, Int_t runIndex, TGraph* graph)
{
    graphArr[cycleIndex]->SetAGraph(runIndex, graph);
}

#endif