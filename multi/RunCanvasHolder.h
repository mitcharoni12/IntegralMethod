#ifndef RUNCANVASHOLDER_H
#define RUNCANVASHOLDER_H

#include <string>
#include "TCanvas.h"

using namespace std;

class RunCanvasHolder{
private:
    TCanvas** canvasArr;
    Int_t lowerRunIndex, numRuns;
public:
    RunCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t cycleIndex, string name);
    ~RunCanvasHolder();
    TCanvas* GetACanvas(Int_t canvasIndex);
};

RunCanvasHolder::RunCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t cycleIndex, string name)
{
    string totalHistoName;
    numRuns = upperRunIndex - lowerRunIndex + 1;

    this->lowerRunIndex = lowerRunIndex;

    //need the plus one becuase if we want run 0 to 0, for our loops and num runs we need the plus one to accuratly calculate
    numRuns = upperRunIndex - lowerRunIndex + 1;
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

//need - lowerRunIndex because if we want to display index 4 through 7, that is a array of size 4
//however, the index of the array is only 0-3, not 4-7, so we have to index a bit weird
//EX: if we want to get canvas 4, that is the 0th element in the array, so for canvas index we get 4 but need 0, so we subtract lower run index
TCanvas* RunCanvasHolder::GetACanvas(Int_t canvasIndex)
{
    canvasIndex = canvasIndex - lowerRunIndex;
    return canvasArr[canvasIndex];
}

#endif