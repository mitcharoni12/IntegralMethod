#ifndef PARAMETERVALUE_H
#define PARAMETERVALUE_H

using namespace std;

/// Holds N0, Half-life, and Decay constant taken from the user for a single element in the decay chain.
/// These values are the inital parameters for the fit functions are the values used to generate the simulated data.
/// Can also take a range for these values if a user wants to set a specific range on the fit values. Currently in the program the fit parameter boundries are just multiples of these values.
class ParameterValue{
    private:
        Double_t lowerRangeHalfLife, upperRangeHalfLife, valueHalfLife;
        bool fixHalfLife;

        Double_t lowerRangeDecayConst, upperRangeDecayConst, valueDecayConst;
        bool fixDecayConst;

        Double_t lowerRangeBatemanN0, upperRangeBatemanN0, batemanN0;
        Double_t lowerRangeIntegralN0, upperRangeIntegralN0, integralN0;
        bool fixBatemanN0, fixIntegralN0;
    public:
        void SetLowerRangeHalfLife(Double_t lower){lowerRangeHalfLife = lower;}
        void SetUpperRangeHalfLife(Double_t upper){upperRangeHalfLife = upper;}
        void SetHalfLife(Double_t val){valueHalfLife = val;}
        void SetFixHalfLife(bool fixHalfLife){this->fixHalfLife = fixHalfLife; this->fixDecayConst = fixHalfLife;}

        void SetLowerRangeDecayConst(Double_t lower){lowerRangeDecayConst = lower;}
        void SetUpperRangeDecayConst(Double_t upper){upperRangeDecayConst = upper;}
        void SetDecayConst(Double_t val){valueDecayConst = val;}

        void SetLowerRangeBatemanN0(Double_t lower){lowerRangeBatemanN0 = lower;}
        void SetUpperRangeBatemanN0(Double_t upper){upperRangeBatemanN0 = upper;}
        void SetBatemanN0(Double_t val){batemanN0 = val;}
        void SetFixBatemanN0(bool fixN0){this->fixBatemanN0 = fixN0;}

        void SetLowerRangeIntegralN0(Double_t lower){lowerRangeIntegralN0 = lower;}
        void SetUpperRangeIntegralN0(Double_t upper){upperRangeIntegralN0 = upper;}
        void SetIntegralN0(Double_t val){integralN0 = val;}
        void SetFixIntegralN0(bool fixN0){this->fixIntegralN0 = fixN0;}

        Double_t GetLowerRangeHalfLife(){return lowerRangeHalfLife;}
        Double_t GetUpperRangeHalfLife(){return upperRangeHalfLife;}
        Double_t GetHalfLife(){return valueHalfLife;}
        Double_t GetLowerRangeHalfLife10Ns(){return lowerRangeHalfLife*1e8;}
        Double_t GetUpperRangeHalfLife10Ns(){return upperRangeHalfLife*1e8;}
        Double_t GetHalfLife10Ns(){return valueHalfLife*1e8;}
        bool GetFixHalfLife(){return fixHalfLife;}

        Double_t GetLowerRangeDecayConst(){return lowerRangeDecayConst;}
        Double_t GetUpperRangeDecayConst(){return upperRangeDecayConst;}
        Double_t GetDecayConst(){return valueDecayConst;}
        Double_t GetLowerRangeDecayConst10Ns(){return lowerRangeDecayConst*1e-8;}
        Double_t GetUpperRangeDecayConst10Ns(){return upperRangeDecayConst*1e-8;}
        Double_t GetDecayConst10Ns(){return valueDecayConst*1e-8;}
        bool GetFixDecayConst(){return fixDecayConst;}

        Double_t GetLowerRangeBatemanN0(){return lowerRangeBatemanN0;}
        Double_t GetUpperRangeBatemanN0(){return upperRangeBatemanN0;}
        Double_t GetBatemanN0(){return batemanN0;}
        bool GetFixBatemanN0(){return fixBatemanN0;}

        Double_t GetLowerRangeIntegralN0(){return lowerRangeIntegralN0;}
        Double_t GetUpperRangeIntegralN0(){return upperRangeIntegralN0;}
        Double_t GetIntegralN0(){return integralN0;}
        bool GetFixIntegralN0(){return fixIntegralN0;}

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
    lowerRangeBatemanN0 = 0;
    upperRangeBatemanN0 = 0;
    batemanN0 = 0;
    lowerRangeIntegralN0 = 0;
    upperRangeIntegralN0 = 0;
    integralN0 = 0;
}
#endif