#ifndef RESULTSTORAGE_H
#define RESULTSTORAGE_H
//class for holding arrays produced by program

using namespace std;

template<class myTemp>
class resultStorage{
    private:
        int arrayDimention, d1Size, d2Size;
        myTemp** doubleArrStorage;
        myTemp* singleArrStorage;
    public:
        resultStorage(int arrayDimention, int d1Size, int d2Size);
        ~resultStorage();
        myTemp** getDoubleArrStorage(){return doubleArrStorage;}
        myTemp* getSingleArrStorage(){return singleArrStorage;}
};

template<typename myTemp>
resultStorage<myTemp>::~resultStorage()
{
    if(arrayDimention == 2)
    {
        for(int i = 0; i < d2Size; i++)
        {
            delete doubleArrStorage[i];
        }
        delete doubleArrStorage;
    }
    if(arrayDimention == 1)
    {
        delete singleArrStorage;
    }
}

template<typename myTemp>
resultStorage<myTemp>::resultStorage(int arrayDimention, int d1Size, int d2Size)
{
    this->arrayDimention = arrayDimention;
    this->d1Size = d1Size;
    this->d2Size = d2Size;
    if(arrayDimention == 2)
    {
        doubleArrStorage = new myTemp* [d2Size];
        for(int i = 0; i < d2Size; i++)
        {
            doubleArrStorage[i] = new myTemp[d1Size];
        }
    }else if(arrayDimention == 1)
    {
        singleArrStorage = new myTemp [d1Size];
    }
}

#endif