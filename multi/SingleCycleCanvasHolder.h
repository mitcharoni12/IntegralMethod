#ifndef SINGLECYCLECANVASHOLDER_H
#define SINGLECYCLECANVASHOLDER_H

#include "TCanvas.h"
#include "CycleCanvasHolder.h"

using namespace std;

class SingleCycleCanvasHolder{
private:
    CycleCanvasHolder** singleCycleCanvasHolder;
    Int_t numElements;
public:
    SingleCycleCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex, Int_t numElements, string* elementNames, string name);
    ~SingleCycleCanvasHolder();
    TCanvas* GetACanvas(Int_t elementIndex, Int_t cycleIndex, Int_t canvasIndex);
};

SingleCycleCanvasHolder::SingleCycleCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t lowerCycleIndex, Int_t upperCycleIndex, Int_t numElements, string* elementNames, string name)
{
    string tempName;

    this->numElements = numElements;
    singleCycleCanvasHolder = new CycleCanvasHolder* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        tempName = elementNames[i] + " " + name;
        singleCycleCanvasHolder[i] = new CycleCanvasHolder(lowerRunIndex, upperRunIndex, lowerCycleIndex, upperCycleIndex, tempName);
    }
}

SingleCycleCanvasHolder::~SingleCycleCanvasHolder()
{
    for(int i = 0; i < numElements; i++)
    {
        delete singleCycleCanvasHolder[i];
    }
    delete [] singleCycleCanvasHolder;
}

TCanvas* SingleCycleCanvasHolder::GetACanvas(Int_t cycleIndex, Int_t canvasIndex, Int_t elementIndex)
{
    return singleCycleCanvasHolder[elementIndex]->GetACanvas(cycleIndex, canvasIndex);
}

#endif