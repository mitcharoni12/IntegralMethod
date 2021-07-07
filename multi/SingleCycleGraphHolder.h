/*
    PURPOSE: used to store the single element integral graphs
*/
#ifndef SINGLECYCLEGRAPHHOLDER_H
#define SINGLECYCLEGRAPHHOLDER_H

#include "TGraph.h"
#include "CycleGraphHolder.h"

using namespace std;

class SingleCycleGraphHolder{
private:
    Int_t numElements;
    CycleGraphHolder** GraphArr;
public:
    SingleCycleGraphHolder(Int_t numCycles, Int_t numElements, Int_t numRuns, string GraphName, Int_t* numPoints, string* elementNames);
    ~SingleCycleGraphHolder();
    TGraph* GetAGraph(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex);
    void SetAGraph(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex, TGraph* Graph);
};

SingleCycleGraphHolder::SingleCycleGraphHolder(Int_t numCycles, Int_t numElements, Int_t numRuns, string GraphName, Int_t* numPoints, string* elementNames)
{
    this->numElements = numElements;
    string GraphParameter = GraphName;
    
    GraphArr = new CycleGraphHolder* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        GraphName = elementNames[i] + " " + GraphName;
        GraphArr[i] = new CycleGraphHolder(numCycles, numRuns, GraphName, numPoints);
        GraphName = GraphParameter;
    }
}

SingleCycleGraphHolder::~SingleCycleGraphHolder()
{
    for(int i = 0; i < numElements; i++)
    {
        delete GraphArr[i];
    }
    delete [] GraphArr;
}

TGraph* SingleCycleGraphHolder::GetAGraph(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex)
{
    return GraphArr[elementIndex]->GetAGraph(cycleIndex, runIndex);
}

void SingleCycleGraphHolder::SetAGraph(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex, TGraph* Graph)
{
    GraphArr[elementIndex]->SetAGraph(cycleIndex, runIndex, Graph);
}

#endif