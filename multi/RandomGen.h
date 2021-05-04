#ifndef RANDOMGEN_H
#define RANDOMGEN_H

#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"

using namespace std;

class RandomGen{
    private:
        Double_t binWidth;
        Double_t numBins;
        Double_t timeRun;
        TRandom3 rand;
        TF1* function;
        //helper function
        void clacBinWidth();
        //setter function
        void setNumBins(Double_t numBins){this->numBins = numBins; clacBinWidth();}
        void setTimeRun(Double_t timeRun){this->timeRun = timeRun; clacBinWidth();}
    public:
        RandomGen(Double_t numBins, Double_t timeRun, TF1* function);
        Double_t genRandomNumber();
};

RandomGen::RandomGen(Double_t numBins, Double_t timeRun, TF1* function)
{
    rand->se
    this->numBins = numBins;
    this->timeRun = timeRun;
    this->function = function;
    clacBinWidth();
}

void RandomGen::clacBinWidth()
{
    this->binWidth = (numBins/timeRun);
}

Double_t RandomGen::genRandomNumber()
{
    Double_t randomHold;
    randomHold = function->GetRandom();
    randomHold = randomHold - (randomHold%binWidth) + (binWidth/2);
    return randomHold;
}

#endif