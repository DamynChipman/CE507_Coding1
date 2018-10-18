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

// Global variables
float DEL_X;

float kElement(int a, int b) {
    Eigen::Matrix2f mat;
    
    float value = 1/DEL_X;
    mat(0,0) = value;
    mat(0,1) = -value;
    mat(1,0) = -value;
    mat(1,1) = value;
    //cout << " In kElement... Returning: " << mat(a,b) << endl;
    return mat(a,b);
}

float fElement(int a) {
    Eigen::Vector2f vec;
    vec(0) = 1;
    vec(1) = 1;
    //cout << " In fElement... Returning: " << vec(a) << endl;
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
    int* ID = (int*)malloc((N + 1) * sizeof(int));
    for (int i = 0; i < (N+1); i++) {
        ID[i] = i;
    }
    ID[N] = -1;
    cout << "ID: " << endl;
    for (int i = 0; i < N+1; i ++) { cout << ID[i] << " ";}
    cout << endl << "IEN: " << endl;
    
    int** IEN = new int*[2];
    for (int i = 0; i < 2; i++) {
        IEN[i] = new int[N];
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0) { IEN[i][j] = ID[j]; }
            else { IEN[i][j] = ID[j]+1; }
            
            cout << IEN[i][j] << " ";
        }
        cout << endl;
    }
    
    cout << "LM: " << endl;
    int** LM = new int*[2];
    for (int i = 0; i < 2; i++) {
        LM[i] = new int[N];
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < N; j++) {
            LM[i][j] = ID[IEN[i][j]];
            
            cout << LM[i][j] << " ";
        }
        cout << endl;
    }

    // Perform FE
    Eigen::Vector3f coefs = FE1D(ID, IEN, LM, kElement, fElement, N);
    cout << endl << coefs;
    
    return 0;
}
