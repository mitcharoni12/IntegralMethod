#ifndef PARAMETERVALUE_H
#define PARAMETERVALUE_H

using namespace std;

/// Holds N0, Half-life, and Decay constant taken from the user for a single element in the decay chain.
/// These values are the inital parameters for the fit functions are the values used to generate the simulated data.
/// Can also take a range for these values if a user wants to set a specific range on the fit values. Currently in the program the fit parameter boundries are just multiples of these values.
class ParameterValue{
    private:
        Double_t lowerRangeHalfLife, upperRangeHalfLife, valueHalfLife;
        bool fixN0;

        Double_t lowerRangeDecayConst, upperRangeDecayConst, valueDecayConst;
        bool fixDecayConst;

        Double_t lowerRangeN0, upperRangeN0, N0;
        bool fixHalfLife;
    public:
        void SetLowerRangeHalfLife(Double_t lower){lowerRangeHalfLife = lower;}
        void SetUpperRangeHalfLife(Double_t upper){upperRangeHalfLife = upper;}
        void SetHalfLife(Double_t val){valueHalfLife = val;}
        void SetFixHalfLife(bool fixHalfLife){this->fixHalfLife = fixHalfLife; this->fixDecayConst = fixHalfLife;}

        void SetLowerRangeDecayConst(Double_t lower){lowerRangeDecayConst = lower;}
        void SetUpperRangeDecayConst(Double_t upper){upperRangeDecayConst = upper;}
        void SetDecayConst(Double_t val){valueDecayConst = val;}

        void SetLowerRangeN0(Double_t lower){lowerRangeN0 = lower;}
        void SetUpperRangeN0(Double_t upper){upperRangeN0 = upper;}
        void SetN0(Double_t val){N0 = val;}
        void SetFixN0(bool fixN0){this->fixN0 = fixN0;}

        Double_t GetLowerRangeHalfLife(){return lowerRangeHalfLife;}
        Double_t GetUpperRangeHalfLife(){return upperRangeHalfLife;}
        Double_t GetHalfLife(){return valueHalfLife;}
        Double_t GetLowerRangeHalfLife10Ns(){return lowerRangeHalfLife*1e8;}
        Double_t GetUpperRangeHalfLife10Ns(){return upperRangeHalfLife*1e8;}
        Double_t GetHalfLife10Ns(){return valueHalfLife*1e8;}

        Double_t GetLowerRangeDecayConst(){return lowerRangeDecayConst;}
        Double_t GetUpperRangeDecayConst(){return upperRangeDecayConst;}
        Double_t GetDecayConst(){return valueDecayConst;}
        Double_t GetLowerRangeDecayConst10Ns(){return lowerRangeDecayConst*1e-8;}
        Double_t GetUpperRangeDecayConst10Ns(){return upperRangeDecayConst*1e-8;}
        Double_t GetDecayConst10Ns(){return valueDecayConst*1e-8;}

        Double_t GetLowerRangeN0(){return lowerRangeN0;}
        Double_t GetUpperRangeN0(){return upperRangeN0;}
        Double_t GetN0(){return N0;}
        ParameterValue();
};

ParameterValue::ParameterValue()
{
    lowerRangeHalfLife = 0;
    upperRangeHalfLife = 0;
    valueHalfLife = 0;
    lowerRangeDecayConst = 0;
    upperRangeDecayConst = 0;
    valueDecayConst = 0;
    lowerRangeN0 = 0;
    upperRangeN0 = 0;
    N0 = 0;
}
#endif