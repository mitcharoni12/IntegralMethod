#ifndef PARAMETERVALUE_H
#define PARAMETERVALUE_H

using namespace std;

/// Holds N0, Half-life, and Decay constant taken from the user for a single element in the decay chain.
/// These values are the inital parameters for the fit functions are the values used to generate the simulated data.
/// Can also take a range for these values if a user wants to set a specific range on the fit values. Currently in the program the fit parameter boundries are just multiples of these values.
class ParameterValue{
    private:
        Double_t lowerRangeHalfLife, upperRangeHalfLife, valueHalfLife;

        Double_t lowerRangeDecayConst, upperRangeDecayConst, valueDecayConst;

        Double_t lowerRangeInitValue, upperRangeInitValue, initValue;
    public:
        void setLowerRangeHalfLife(Double_t lower){lowerRangeHalfLife = lower;}
        void setUpperRangeHalfLife(Double_t upper){upperRangeHalfLife = upper;}
        void setValueHalfLife(Double_t val){valueHalfLife = val;}

        void setLowerRangeDecayConst(Double_t lower){lowerRangeDecayConst = lower;}
        void setUpperRangeDecayConst(Double_t upper){upperRangeDecayConst = upper;}
        void setValueDecayConst(Double_t val){valueDecayConst = val;}

        void setLowerRangeInitValue(Double_t lower){lowerRangeInitValue = lower;}
        void setUpperRangeInitValue(Double_t upper){upperRangeInitValue = upper;}
        void setInitValue(Double_t val){initValue = val;}

        Double_t getLowerRangeHalfLife(){return lowerRangeHalfLife;}
        Double_t getUpperRangeHalfLife(){return upperRangeHalfLife;}
        Double_t getValueHalfLife(){return valueHalfLife;}
        Double_t getLowerRangeHalfLife10Ns(){return lowerRangeHalfLife*1e8;}
        Double_t getUpperRangeHalfLife10Ns(){return upperRangeHalfLife*1e8;}
        Double_t getValueHalfLife10Ns(){return valueHalfLife*1e8;}

        Double_t getLowerRangeDecayConst(){return lowerRangeDecayConst;}
        Double_t getUpperRangeDecayConst(){return upperRangeDecayConst;}
        Double_t getValueDecayConst(){return valueDecayConst;}
        Double_t getLowerRangeDecayConst10Ns(){return lowerRangeDecayConst*1e-8;}
        Double_t getUpperRangeDecayConst10Ns(){return upperRangeDecayConst*1e-8;}
        Double_t getValueDecayConst10Ns(){return valueDecayConst*1e-8;}

        Double_t getLowerRangeInitValue(){return lowerRangeInitValue;}
        Double_t getUpperRangeInitValue(){return upperRangeInitValue;}
        Double_t getInitValue(){return initValue;}
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
    lowerRangeInitValue = 0;
    upperRangeInitValue = 0;
    initValue = 0;
}
#endif