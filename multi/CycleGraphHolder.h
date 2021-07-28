#ifndef CYCLEGRAPHHOLDER_H
#define CYCLEGRAPHHOLDER_H

using namespace std;

#include "TGraph.h"
#include "RunGraphHolder.h"

/// Used to store the graphs for the multiple cycles of a single element. These graphs are derived from a coresponding integral histogram.
/// EX: If we are doing 20 cycles of 20 runs for the decay chain 144Cs->144Ba->144La this class would hold all the histograms for the 20 cycles of La, each cycle held by a RunHistoHolder class.
class CycleGraphHolder{
private:
    Int_t numCycles;
    RunGraphHolder** graphArr; ///< Holds all the graphss, array of RunGraphHolder objects.
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

/// Gets a graph at a specific cycle and run index.
TGraph* CycleGraphHolder::GetAGraph(Int_t cycleIndex, Int_t runIndex)
{
    return graphArr[cycleIndex]->GetAGraph(runIndex);
}

/// Sets a graph at a specific cycle and run index.
void CycleGraphHolder::SetAGraph(Int_t cycleIndex, Int_t runIndex, TGraph* graph)
{
    graphArr[cycleIndex]->SetAGraph(runIndex, graph);
}

#endif