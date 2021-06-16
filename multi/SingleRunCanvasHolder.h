#ifndef SINGLERUNCANVASHOLDER_H
#define SINGLERUNCANVASHOLDER_H

#include "TCanvas.h"
#include "RunCanvasHolder.h"

using namespace std;

class SingleRunCanvasHolder{
private:
    RunCanvasHolder** singleCanvasArr;
    Int_t numElements;
public:
    SingleRunCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t cycleIndex, string name, Int_t numElements, string* elementNames);
    ~SingleRunCanvasHolder();
    TCanvas* GetACanvas(Int_t elementIndex, Int_t canvasIndex);
};

SingleRunCanvasHolder::SingleRunCanvasHolder(Int_t lowerRunIndex, Int_t upperRunIndex, Int_t cycleIndex, string name, Int_t numElements, string* elementNames)
{
    string tempName;
    this->numElements = numElements;

    singleCanvasArr = new RunCanvasHolder* [numElements];
    for(int i = 0; i < numElements; i++)
    {
        tempName = elementNames[i] + " " + name;
        singleCanvasArr[i] = new RunCanvasHolder(lowerRunIndex, upperRunIndex, cycleIndex, tempName); 
    }
}

SingleRunCanvasHolder::~SingleRunCanvasHolder()
{
    for(int i = 0; i < numElements; i++)
    {
        delete singleCanvasArr[i];
    }
    delete [] singleCanvasArr;
}

TCanvas* SingleRunCanvasHolder::GetACanvas(Int_t elementIndex, Int_t canvasIndex)
{
    return singleCanvasArr[elementIndex]->GetACanvas(canvasIndex);
}

#endif