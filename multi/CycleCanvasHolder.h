#ifndef CYCLECANVASHOLDER_H
#define CYCLECANVASHOLDER_H

#include "RunCanvasHolder.h"
#include "TCanvas.h"

using namespace std;

class CycleCanvasHolder{
private:
    RunCanvasHolder** cycleCanvasHolder;
    Int_t numCycles, lowerCycleIndex;
public:
    CycleCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex, string name);
    ~CycleCanvasHolder();
    TCanvas* GetACanvas(Int_t cycleIndex, Int_t canvasIndex);
};

CycleCanvasHolder::CycleCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex, string name)
{   //have to add one because if we want just cycle 0, that would be 0 to 0 and for our looks and num cycles we need the plus one
    numCycles = upperCycleIndex - lowerCycleIndex + 1;
    this->lowerCycleIndex = lowerCycleIndex;

    cycleCanvasHolder = new RunCanvasHolder* [numCycles];
    for(int i = 0; i < numCycles; i++)
    {
        cycleCanvasHolder[i] = new RunCanvasHolder(lowerRunIndex, upperRunIndex, i+1, name);
    }
}

CycleCanvasHolder::~CycleCanvasHolder()
{
    for(int i = 0; i < numCycles; i++)
    {
        delete cycleCanvasHolder[i];
    }
    delete [] cycleCanvasHolder;
}

//need - loweCycleIndex because if we want to display index 4 through 7, that is a array of size 4
//however, the index of the array is only 0-3, not 4-7, so we have to index a bit weird
//EX: if we want to get object 4, that is the 0th element in the array, so for cycle index we get 4 but need 0, so we subtract lower run index
TCanvas* CycleCanvasHolder::GetACanvas(Int_t cycleIndex, Int_t canvasIndex)
{
    cycleIndex = cycleIndex - lowerCycleIndex;
    return cycleCanvasHolder[cycleIndex]->GetACanvas(canvasIndex);
}

#endif