#include <iomanip>
#include <utility>
#include <vector>
#include <cstdio>
#include <sys/stat.h>
#include "TGraph.h"
#include "TDirectory.h"

using namespace std;

void DisplayFile()
{
    TFile* in_file = new TFile("outputFile.root");
    TDriectory
    TGraph **graphArr = new TGraph* [12];
    TGraph* graphBoi;
    for(int i = 0; i < 1; i++)
    {
        in_file->GetObject((("Graph;" + to_string(i+1)).c_str()), graphArr[i]);
    }
    
    TCanvas *regularCSErrorCanvas = new TCanvas("regularCSErrorCanvas", "regularCSErrorCanvas", 500, 500);
    /*
    TCanvas *integralCSErrorCanvas = new TCanvas("integralCSErrorCanvas", "integralCSErrorCanvas", 500, 500);
    TCanvas *regularCSParaCanvas = new TCanvas("regularCSParaCanvas", "regularCSParaCanvas", 500, 500);
    TCanvas *integralCSParaCanvas = new TCanvas("integralCSParaCanvas", "integralCSParaCanvas", 500, 500);
    TCanvas *CS0regularValueCanvas = new TCanvas("CS0regularValueCanvas", "CS0regularValueCanvas", 500, 500);
    TCanvas *CS0integralValueCanvas = new TCanvas("CS0integralValueCanvas", "CS0integralValueCanvas", 500, 500);
    TCanvas *regularBAErrorCanvas = new TCanvas("regularBAErrorCanvas", "regularBAErrorCanvas", 500, 500);
    TCanvas *integralBAErrorCanvas = new TCanvas("integralBAErrorCanvas", "integralBAErrorCanvas", 500, 500);
    TCanvas *regularBAParaCanvas = new TCanvas("regularBAParaCanvas", "regularBAParaCanvas", 500, 500);
    TCanvas *integralBAParaCanvas = new TCanvas("integralBAParaCanvas", "integralBAParaCanvas", 500, 500);
    TCanvas *BA0regularValueCanvas = new TCanvas("BA0regularValueCanvas", "BA0regularValueCanvas", 500, 500);
    TCanvas *BA0integralValueCanvas = new TCanvas("BA0integralValueCanvas", "BA0integralValueCanvas", 500, 500);
    */
    //regularCSErrorCanvas->cd();
    //graphArr[0]->Draw();
    in_file->Close();
}