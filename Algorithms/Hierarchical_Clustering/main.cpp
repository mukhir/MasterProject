#include <iostream>
#include <fstream>

#include "HierarchicalClustering.hpp"



int main()
{   
    int outputCount = 21000;
    splatPoint* sArr = new splatPoint[outputCount];
    
    HierClustering hClust;
    hClust.init("armadillo.points");
	hClust.contract(sArr, outputCount);
    
    std::ofstream wfile("output");
    wfile.write((char*)sArr, outputCount*sizeof(splatPoint));
    wfile.close();
    
    delete[] sArr;
    
    return 0;
}