#ifndef PARAMETERVALUE_H
#define PARAMETERVALUE_H

using namespace std;

/// Holds N0, Half-life, and Decay constant taken from the user for a single element in the decay chain.
/// These values are the inital parameters for the fit functions are the values used to generate the simulated data.
/// Can also take a range for these values if a user wants to set a specific range on the fit values. Currently in the program the fit parameter boundries are just multiples of these values.
class ParameterValue{
    private:
        Double_t lowerRangeHalfLife;
        Double_t upperRangeHalfLife;
        Double_t rangeAverageHalfLife;
        Double_t valueHalfLife;
        bool isValueHalfLife; ///< True = Half-life is a direct value rather than a value range.

        Double_t lowerRangeDecayConst;
        Double_t upperRangeDecayConst;
        Double_t rangeAverageDecayConst;
        Double_t valueDecayConst;
        bool isValueDecayConst; ///< True = Decay constant is a direct value rather than a value range.

        Double_t lowerRangeInitValue;
        Double_t upperRangeInitValue;
        Double_t rangeAverageInitValue;
        Double_t initValue;
        bool isValueInitValue; ///< True = Inital value is a direct value rather than a value range.
    public:
        void setLowerRangeHalfLife(Double_t lower){lowerRangeHalfLife = lower;}
        void setUpperRangeHalfLife(Double_t upper){upperRangeHalfLife = upper;}
        void setValueHalfLife(Double_t val){valueHalfLife = val;}
        void setIsValueHalfLife(bool isVal){isValueHalfLife = isVal;}

        void setLowerRangeDecayConst(Double_t lower){lowerRangeDecayConst = lower;}
        void setUpperRangeDecayConst(Double_t upper){upperRangeDecayConst = upper;}
        void setValueDecayConst(Double_t val){valueDecayConst = val;}
        void setIsValueDecayConst(bool isVal){isValueDecayConst = isVal;}

        void setLowerRangeInitValue(Double_t lower){lowerRangeInitValue = lower;}
        void setUpperRangeInitValue(Double_t upper){upperRangeInitValue = upper;}
        void setInitValue(Double_t val){initValue = val;}
        void setIsValueInitValue(bool isVal){isValueInitValue = isVal;}

        Double_t getLowerRangeHalfLife(){return lowerRangeHalfLife;}
        Double_t getUpperRangeHalfLife(){return upperRangeHalfLife;}
        Double_t getRangeAverageHalfLife(){return rangeAverageHalfLife;}
        Double_t getValueHalfLife(){return valueHalfLife;}
        bool getIsValueHalfLife(){return isValueHalfLife;}

        Double_t getLowerRangeDecayConst(){return lowerRangeDecayConst;}
        Double_t getUpperRangeDecayConst(){return upperRangeDecayConst;}
        Double_t getRangeAverageDecayConst(){return rangeAverageDecayConst;}
        Double_t getValueDecayConst(){return valueDecayConst;}
        bool getIsValueDecayConst(){return isValueDecayConst;}

        Double_t getLowerRangeInitValue(){return lowerRangeInitValue;}
        Double_t getUpperRangeInitValue(){return upperRangeInitValue;}
        Double_t getRangeAverageInitValue(){return rangeAverageInitValue;}
        Double_t getInitValue(){return initValue;}
        bool getIsValueInitValue(){return isValueInitValue;}
        ParameterValue();
        void calculateRangeAverageHalfLife();
        void calculateRangeAverageDecayConst();
        void calculateRangeAverageInitValue();
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

void ParameterValue::calculateRangeAverageHalfLife()
{
    rangeAverageHalfLife = ((lowerRangeHalfLife + upperRangeHalfLife) / 2);
}

void ParameterValue::calculateRangeAverageDecayConst()
{
    rangeAverageDecayConst = ((lowerRangeDecayConst + upperRangeDecayConst) / 2);
}

void ParameterValue::calculateRangeAverageInitValue()
{
    rangeAverageInitValue = ((lowerRangeInitValue + upperRangeInitValue) / 2);
}

#endif