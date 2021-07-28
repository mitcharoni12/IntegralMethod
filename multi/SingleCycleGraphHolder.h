#ifndef SINGLECYCLEGRAPHHOLDER_H
#define SINGLECYCLEGRAPHHOLDER_H

#include "TGraph.h"
#include "CycleGraphHolder.h"

using namespace std;

/// Holds all the graphs for all the single histograms of the program.
/// EX: If we are doing 20 cycles of 20 runs for the decay chain 144Cs->144Ba->144La we have single graphs for each element in the decay chain and a total integral graph.
/// The total integral graph for the 20 cycles of 20 runs are held by the CylceGraphHolder class but the single graphs for Cs, Ba, and La are all held by this class.
class SingleCycleGraphHolder{
private:
    Int_t numElements;
    CycleGraphHolder** GraphArr; ///< array of size numElements of CycleHistoHolder objects.
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

/// Gets a graph for a specific run, cycle, and element index.
TGraph* SingleCycleGraphHolder::GetAGraph(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex)
{
    return GraphArr[elementIndex]->GetAGraph(cycleIndex, runIndex);
}

/// Sets a graph for a specific run, cycle, and element index.
void SingleCycleGraphHolder::SetAGraph(Int_t cycleIndex, Int_t runIndex, Int_t elementIndex, TGraph* Graph)
{
    GraphArr[elementIndex]->SetAGraph(cycleIndex, runIndex, Graph);
}

#endif