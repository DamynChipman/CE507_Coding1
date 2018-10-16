//   Created by Damyn Chipman on 10/6/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   main.cpp
//   PROJECT:   CE507_Coding1

#include <iostream>
#include <Eigen/Dense>
//#include "/Users/Damyn/Documents/Coding Files/XCode/Numerical Data Structures/Numerical Data Structures/Matrix.h"
#include "Linear.h"
#include "FE1D.h"

using namespace std;
using namespace Eigen;

float DEL_X;

float kElement(int a, int b) {
    Eigen::Matrix3f mat;

    float value = 1/DEL_X;
    mat(0,0) = value;
    mat(0,1) = -value;
    mat(1,0) = -value;
    mat(1,1) = value;
    
    return mat(a,b);
}

float fElement(int a) {
    Eigen::Vector3f vec;
    // TODO: Add f_e vector
    return vec(a);
}

int main(int argc, const char * argv[]) {
    
    // Create function and plotting domain
    float xL = 0;
    float xR = 1;
    int N = 10;
    int NPlot = 10*N;
    Linear Nfunction(xL,xR,N);
    Linear::Domain domain = Nfunction.getDomain();
    Linear::Domain plotDomain(xL,xR,NPlot);
    DEL_X = domain.del_x_;
    
    // Create array mappings
    int* ID;
    int** IEN;
    int** LM;
    
    // Create element matrices (functions)

    // Perform FE
    Eigen::Vector3f coefs = FE1D(ID, IEN, LM, kElement, fElement, N);
    
    return 0;
}
