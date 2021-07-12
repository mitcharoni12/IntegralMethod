/*
    PURPOSE: used to store the graphs for the integral histograms
*/
#ifndef RUNGRAPHHOLDER_H
#define RUNGRAPHHOLDER_H

#include "TGraph.h"

using namespace std;

class RunGraphHolder{
private:
    TGraph** graphArr;
public:
    RunGraphHolder(Int_t cycleNum, Int_t numRuns, string GraphName, Int_t numPoints);
    ~RunGraphHolder();
    TGraph* GetAGraph(Int_t GraphIndex);
    void SetAGraph(Int_t runIndex, TGraph* Graph);
};

RunGraphHolder::RunGraphHolder(Int_t cycleNum, Int_t numRuns, string graphName, Int_t numPoints)
{
    string tempGraphName;

    graphArr = new TGraph* [numRuns];
    for(int i = 0; i < numRuns; i++)
    {
        cout << "RUN " << i << " " << tempGraphName << " Cycle: " << cycleNum << " Runs: " << numRuns << " points " << numPoints << endl;
        tempGraphName = graphName + " Run: " + (to_string(i+1)).c_str() + " Cycle: " + (to_string(cycleNum+1)).c_str();
        graphArr[i] = new TGraph(numPoints);
        graphArr[i]->SetNameTitle(tempGraphName.c_str());
    }
}

RunGraphHolder::~RunGraphHolder()
{
    delete [] graphArr;
}

TGraph* RunGraphHolder::GetAGraph(Int_t runIndex)
{
    TGraph* tempGraph = graphArr[runIndex];
    return tempGraph;
}

void RunGraphHolder::SetAGraph(Int_t runIndex, TGraph* Graph)
{
    graphArr[runIndex] = Graph;
}

#endif