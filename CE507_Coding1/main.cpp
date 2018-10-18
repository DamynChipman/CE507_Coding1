//   Created by Damyn Chipman on 10/6/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   main.cpp
//   PROJECT:   CE507_Coding1

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "Linear.h"
#include "FE1D.h"

using namespace std;
using namespace Eigen;

// Global variables
float DEL_X; // To be set later with specific discrete domain

float kElement(int a, int b) {
    Eigen::Matrix2f mat;
    
    float value = 1/DEL_X;
    mat(0,0) = value;
    mat(0,1) = -value;
    mat(1,0) = -value;
    mat(1,1) = value;
    
    return mat(a,b);
}

float fElement1(int a) {
    Eigen::Vector2f vec;
    float c = 1;
    float f1 = c;
    float f2 = c;
    
    float value = DEL_X/6;
    vec(0) = value*(2*f1 + f2);
    vec(1) = value*(f1 + 2*f2);
    
    return vec(a);
}

float fElement2(int a) {
    Eigen::Vector2f vec;
    
    float f1 = 0;
    float f2 = 1;
    // TODO: Fix these values
    float value = DEL_X/6;
    vec(0) = value*(2*f1 + f2);
    vec(1) = value*(f1 + 2*f2);
    
    return vec(a);
}

float fElement3(int a) {
    Eigen::Vector2f vec;
    // TODO: Fix these values
    float f1 = 0;
    float f2 = 1;
    
    float value = DEL_X/6;
    vec(0) = value*(2*f1 + f2);
    vec(1) = value*(f1 + 2*f2);
    
    return vec(a);
}

int main(int argc, const char * argv[]) {
    
    // Output file names and N nodes arrays
    string xNames[4] = {"xN1.csv","xN2.csv","xN3.csv","xN4.csv"};
    string u1Names[4] = {"u1N1.csv","u1N2.csv","u1N3.csv","u1N4.csv"};
    string u2Names[4] = {"u2N1.csv","u2N2.csv","u2N3.csv","u2N4.csv"};
    string u3Names[4] = {"u3N1.csv","u3N2.csv","u3N3.csv","u3N4.csv"};
    int N[4] = {10,100,1000,10000};
    
    // --- BEGIN LOOP ---
    cout << " --- BEGINNING PROGRAM LOOP --- " << endl;
    for (int n = 0; n < 1; n++) {
        cout << " n = " << n << " ---------- N = " << N[n] << " ------------ " << endl;
        // Create function and plotting domain
        cout << "   generating domains... " << endl;
        float xL = 0;
        float xR = 1;
        int NPlot = 100*N[n];
        Linear Nfunction(xL,xR,N[n]);
        Linear::Domain domain = Nfunction.getDomain();
        Linear::Domain plotDomain(xL,xR,NPlot);
        DEL_X = domain.del_x_;
        
        // Create array mappings
        cout << "   generating array mappings... " << endl;
        int* ID = (int*)malloc((N[n] + 1) * sizeof(int));
        for (int i = 0; i < (N[n]+1); i++) {
            ID[i] = i;
        }
        ID[N[n]] = -1;
        
        int** IEN = new int*[2];
        for (int i = 0; i < 2; i++) {
            IEN[i] = new int[N[n]];
        }
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < N[n]; j++) {
                if (i == 0) { IEN[i][j] = ID[j]; }
                else { IEN[i][j] = ID[j]+1; }
            }
        }
        
        int** LM = new int*[2];
        for (int i = 0; i < 2; i++) {
            LM[i] = new int[N[n]];
        }
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < N[n]; j++) {
                LM[i][j] = ID[IEN[i][j]];
            }
        }

        // Perform FE
        cout << "   calculating coefs1... " << endl;
        Eigen::VectorXf coefs1 = FE1D(ID, IEN, LM, kElement, fElement1, N[n]); // For f = c
        cout << "   calculating coefs2... " << endl;
        Eigen::VectorXf coefs2 = FE1D(ID, IEN, LM, kElement, fElement2, N[n]); // For f = x
        cout << "   calculating coefs3... " << endl;
        Eigen::VectorXf coefs3 = FE1D(ID, IEN, LM, kElement, fElement3, N[n]); // For f = x^2
        
        // Generate solution
        cout << "   generating solutions... " << endl;
        float* u_h1 = (float*)malloc(NPlot * sizeof(float));
        float* u_h2 = (float*)malloc(NPlot * sizeof(float));
        float* u_h3 = (float*)malloc(NPlot * sizeof(float));
        for (int i = 0; i < NPlot; i++) {
            u_h1[i] = 0;
            u_h2[i] = 0;
            u_h3[i] = 0;
            for (int j = 0; j < N[n]; j++) {
                u_h1[i] = u_h1[i] + coefs1(j)*Nfunction.eval(j, plotDomain.x_h_[i]);
                u_h2[i] = u_h2[i] + coefs2(j)*Nfunction.eval(j, plotDomain.x_h_[i]);
                u_h3[i] = u_h3[i] + coefs3(j)*Nfunction.eval(j, plotDomain.x_h_[i]);
            }
        }
        
        // Output results
        cout << "   outputing results... " << endl;
        ofstream x(xNames[n]);
        ofstream u1(u1Names[n]);
        ofstream u2(u2Names[n]);
        ofstream u3(u3Names[n]);
        char delim = ',';
        
        for (int i = 0; i < NPlot; i++) {
            x << plotDomain.x_h_[i] << delim;
            u1 << u_h1[i] << delim;
            u2 << u_h2[i] << delim;
            u2 << u_h3[i] << delim;
        }
        
        x.close();
        u1.close();
        u2.close();
        u3.close();
    }
    
    return 0;
}
