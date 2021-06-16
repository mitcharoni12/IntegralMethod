#ifndef RUNCANVASHOLDER_H
#define RUNCANVASHOLDER_H

#include <string>
#include "TCanvas.h"

using namespace std;

class RunCanvasHolder{
private:
    TCanvas** canvasArr;
public:
    RunCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t cycleIndex, string name);
    ~RunCanvasHolder();
    TCanvas* GetACanvas(Int_t canvasIndex);
};

RunCanvasHolder::RunCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t cycleIndex, string name)
{
    string totalHistoName;
    Int_t numRuns;

    numRuns = upperRunIndex - lowerRunIndex;
    canvasArr = new TCanvas* [numRuns];
    for(int i = 0; i < numRuns; i++)
    {
        totalHistoName = name + " Run: " + (to_string(i + 1 + lowerRunIndex)).c_str() + " Cycle: " + (to_string(cycleIndex)).c_str();
        canvasArr[i] = new TCanvas(totalHistoName.c_str(), totalHistoName.c_str(), 500, 500);
    }
}

RunCanvasHolder::~RunCanvasHolder()
{
    delete [] canvasArr;
}

TCanvas* RunCanvasHolder::GetACanvas(Int_t canvasIndex)
{
    return canvasArr[canvasIndex];
}

#endif