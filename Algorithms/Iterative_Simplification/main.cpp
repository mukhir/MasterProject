#include <iostream>
#include <fstream>

#include "pauly.hpp"



int main()
{   
    int outputCount = 32000;
    splatPoint* sArr = new splatPoint[outputCount];
    
    PaulySimp paulyS;
    paulyS.init("armadillo.points");
	paulyS.contract(sArr, outputCount);
    
    std::ofstream wfile("output");
    wfile.write((char*)sArr, outputCount*sizeof(splatPoint));
    wfile.close();
    
    delete[] sArr;
    
    return 0;
}