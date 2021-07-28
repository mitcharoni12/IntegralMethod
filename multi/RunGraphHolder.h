#ifndef RUNGRAPHHOLDER_H
#define RUNGRAPHHOLDER_H

#include "TGraph.h"

using namespace std;

/// Used to store all the graphs for the multiple runs in a single cycle for a single element. These graphs are derived from a coresponding integral histogram.
/// EX: If we are doing 20 cycles and 20 runs per cycles of the decay chain 144Cs->144Ba->144La. This would hold the 20 graphs of a single cycle for the element La.
class RunGraphHolder{
private:
    TGraph** graphArr;
public:
    RunGraphHolder(Int_t cycleNum, Int_t numRuns, string GraphName, Int_t numPoints);
    ~RunGraphHolder();
    TGraph* GetAGraph(Int_t GraphIndex);
    void SetAGraph(Int_t runIndex, TGraph* Graph);
};

/// Constructor for class. Dynaically allocates graphs and names them properly.
RunGraphHolder::RunGraphHolder(Int_t cycleNum, Int_t numRuns, string graphName, Int_t numPoints)
{
    string tempGraphName;

    graphArr = new TGraph* [numRuns];
    for(int i = 0; i < numRuns; i++)
    {
        tempGraphName = graphName + " Run: " + (to_string(i+1)).c_str() + " Cycle: " + (to_string(cycleNum+1)).c_str();
        graphArr[i] = new TGraph(numPoints);
        graphArr[i]->SetNameTitle(tempGraphName.c_str());
    }
}

RunGraphHolder::~RunGraphHolder()
{
    delete [] graphArr;
}

/// Gets a graph at a specific run index.
TGraph* RunGraphHolder::GetAGraph(Int_t runIndex)
{
    TGraph* tempGraph = graphArr[runIndex];
    return tempGraph;
}

/// Sets a graph at a specific run index.
void RunGraphHolder::SetAGraph(Int_t runIndex, TGraph* Graph)
{
    graphArr[runIndex] = Graph;
}

#endif