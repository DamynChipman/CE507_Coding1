//   Created by Damyn Chipman on 10/6/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   main.cpp
//   PROJECT:   CE507_Coding1

#include <iostream>
#include "/Users/Damyn/Documents/Coding Files/XCode/Numerical Data Structures/Numerical Data Structures/Matrix.h"
#include "Linear.h"
#include "FE1D.h"

using namespace std;



int main(int argc, const char * argv[]) {
    
    // Create function and plotting domain
    float xL = 0;
    float xR = 1;
    int N = 10;
    int NPlot = 10*N;
    Linear Nfunction(xL,xR,N);
    Linear::Domain domain = Nfunction.getDomain();
    Linear::Domain plotDomain(xL,xR,NPlot);
    
    // Create array mappings
    int* ID = (int*)malloc((N-1) * sizeof(int));
    
    
    return 0;
}
