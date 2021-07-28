#ifndef SINGLECYCLECANVASHOLDER_H
#define SINGLECYCLECANVASHOLDER_H

#include "TCanvas.h"
#include "CycleCanvasHolder.h"

using namespace std;

/// Used to store the canvases which will display each individual histogram if the user wishes to display each individual histogram.
/// This class is used to store the canvases which will display the single element histograms.
/// EX: If we are doing 20 cycles of 20 runs of the decay chain 144Cs->144Ba->144La we will have single element histograms and total and integral and Bateman histograms to display.
/// The total integral and bateman canvases for the 20 cycles and runs are stored by the CycleCanvasHolder class. This class holds single element canvases for each element in the decay chain.
class SingleCycleCanvasHolder{
private:
    CycleCanvasHolder** singleCycleCanvasHolder; ///< array of size numElements of CycleCanvasHolder objects.
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

/// Get a canvas at a specific cycle, run and element index.
TCanvas* SingleCycleCanvasHolder::GetACanvas(Int_t cycleIndex, Int_t canvasIndex, Int_t elementIndex)
{
    return singleCycleCanvasHolder[elementIndex]->GetACanvas(cycleIndex, canvasIndex);
}

#endif