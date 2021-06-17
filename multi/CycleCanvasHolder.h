#ifndef CYCLECANVASHOLDER_H
#define CYCLECANVASHOLDER_H

#include "RunCanvasHolder.h"
#include "TCanvas.h"

using namespace std;

class CycleCanvasHolder{
private:
    RunCanvasHolder** cycleCanvasHolder;
    Int_t numCycles;
public:
    CycleCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex, string name);
    ~CycleCanvasHolder();
    TCanvas* GetACanvas(Int_t cycleIndex, Int_t canvasIndex);
};

CycleCanvasHolder::CycleCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex, string name)
{   //have to add one because if we want just cycle 0, that would be 0 to 0 and for our looks and num cycles we need the plus one
    numCycles = upperCycleIndex - lowerCycleIndex + 1;
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

TCanvas* CycleCanvasHolder::GetACanvas(Int_t cycleIndex, Int_t canvasIndex)
{
    cout << "IN GET A CANVAS cycle index: " << cycleIndex << " run index: " << canvasIndex << endl;
    return cycleCanvasHolder[cycleIndex]->GetACanvas(canvasIndex);
}

#endif