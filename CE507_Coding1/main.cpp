//   Created by Damyn Chipman on 10/6/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   main.cpp
//   PROJECT:   CE507_Coding1

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <Eigen/Dense>
#include "functions.h"
#include "Linear.h"
#include "Domain.h"
#include "FE1D.h"

using namespace std;
using namespace Eigen;

// Global variables
int N[4] = {10,100,1000,10000}; // Number of nodes
string XNAMES[4] = {"xN1.csv","xN2.csv","xN3.csv","xN4.csv"};
string U1NAMES[4] = {"u1N1.csv","u1N2.csv","u1N3.csv","u1N4.csv"};
string U2NAMES[4] = {"u2N1.csv","u2N2.csv","u2N3.csv","u2N4.csv"};
string U3NAMES[4] = {"u3N1.csv","u3N2.csv","u3N3.csv","u3N4.csv"};
string UACT1NAMES[4] = {"uAct1N1.csv","uAct1N2.csv","uAct1N3.csv","uAct1N4.csv"};
string UACT2NAMES[4] = {"uAct2N1.csv","uAct2N2.csv","uAct2N3.csv","uAct2N4.csv"};
string UACT3NAMES[4] = {"uAct3N1.csv","uAct3N2.csv","uAct3N3.csv","uAct3N4.csv"};
string IDNAMES[4] = {"IDN1.csv","IDN2.csv","IDN3.csv","IDN4.csv"};
string IENNAMES[4] = {"IENN1.csv","IENN2.csv","IENN3.csv","IENN4.csv"};
string LMNAMES[4] = {"LMN1.csv","LMN2.csv","LMN3.csv","LMN4.csv"};

// ----- Main Function -----
/**
 * @function main
 * @brief Coding Assignment 1 - CE507
 * @returns 1 : int : Program ran successfully
 */
int main() {
    
    // Domain and BC
    float xL = 0;
    float xR = 1;
    float g = 0; // Necessary BC @ x = xR = 1.0
    
    // --- BEGIN LOOP ---
    cout << " --- BEGINNING PROGRAM LOOP --- " << endl;
    clock_t t;   // Timing
    t = clock(); // Timing
    for (int n = 0; n < 4; n++) {
        cout << " n = " << n << " ---------- N = " << N[n] << " ------------ " << endl;
        
        // Create function and plotting domain
        cout << "   generating domains... " << endl;
        int NPlot = N[n]+1;
        Linear Nfunction(xL,xR,N[n]);
        Linear::Domain domain = Nfunction.getDomain();
        Linear::Domain plotDomain(xL,xR,NPlot);
        DEL_X = domain.getDelX();
        
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
        
        // Perform FE1D
        cout << "   calculating coefs1... " << endl;
        Eigen::VectorXf coefs1 = FE1D(ID, IEN, LM, kElement, fElement1, N[n]); // For f = c
        cout << "   calculating coefs2... " << endl;
        Eigen::VectorXf coefs2 = FE1D(ID, IEN, LM, kElement, fElement2, N[n]); // For f = x
        cout << "   calculating coefs3... " << endl;
        Eigen::VectorXf coefs3 = FE1D(ID, IEN, LM, kElement, fElement3, N[n]); // For f = x^2
        
        // Generate solution
        cout << "   generating solutions... " << endl;
        float u_h1[NPlot];
        float u_h2[NPlot];
        float u_h3[NPlot];
        float u_act1[NPlot];
        float u_act2[NPlot];
        float u_act3[NPlot];
        for (int i = 0; i < NPlot; i++) {
            u_h1[i] = 0;
            u_h2[i] = 0;
            u_h3[i] = 0;
            u_act1[i] = uActual(plotDomain.getPoints()[i], 1);
            u_act2[i] = uActual(plotDomain.getPoints()[i], 2);
            u_act3[i] = uActual(plotDomain.getPoints()[i], 3);
            for (int j = 0; j < N[n]; j++) {
                u_h1[i] = u_h1[i] + coefs1(j)*Nfunction.eval(j, plotDomain.getPoints()[i]);
                u_h2[i] = u_h2[i] + coefs2(j)*Nfunction.eval(j, plotDomain.getPoints()[i]);
                u_h3[i] = u_h3[i] + coefs3(j)*Nfunction.eval(j, plotDomain.getPoints()[i]);
            }
        }
        // Add BC (known points)
        u_h1[NPlot] = g;
        u_h2[NPlot] = g;
        u_h3[NPlot] = g;
        u_act1[NPlot-1] = g;
        u_act2[NPlot-1] = g;
        u_act3[NPlot-1] = g;
        
        // Output results
        cout << "   outputing results... " << endl;
        ofstream x(XNAMES[n]);
        ofstream u1(U1NAMES[n]);
        ofstream u2(U2NAMES[n]);
        ofstream u3(U3NAMES[n]);
        ofstream uAct1(UACT1NAMES[n]);
        ofstream uAct2(UACT2NAMES[n]);
        ofstream uAct3(UACT3NAMES[n]);
        ofstream IDOutput(IDNAMES[n]);
        ofstream IENOutput(IENNAMES[n]);
        ofstream LMOutput(LMNAMES[n]);
        char delim = ',';
        
        for (int i = 0; i < NPlot-1; i++) {
            x << plotDomain.getPoints()[i] << delim;
            u1 << u_h1[i] << delim;
            u2 << u_h2[i] << delim;
            u3 << u_h3[i] << delim;
            uAct1 << u_act1[i] << delim;
            uAct2 << u_act2[i] << delim;
            uAct3 << u_act3[i] << delim;
        }
        x << xR << delim;
        u1 << u_h1[NPlot] << delim;
        u2 << u_h2[NPlot] << delim;
        u3 << u_h3[NPlot] << delim;
        uAct1 << u_act1[NPlot] << delim;
        uAct2 << u_act2[NPlot] << delim;
        uAct3 << u_act3[NPlot] << delim;
        
        for (int i = 0; i < N[n]; i++) {
            IDOutput << ID[i] << delim;
            for (int j = 0; j < 2; j++) {
                IENOutput << IEN[j][i] << delim;
                LMOutput << LM[j][i] << delim;
            }
            IENOutput << endl;
            LMOutput << endl;
        }
        
        x.close();
        u1.close();
        u2.close();
        u3.close();
        uAct1.close();
        uAct2.close();
        uAct3.close();
        IDOutput.close();
        IENOutput.close();
        LMOutput.close();
        
        // Memory
        free(ID);
        for (int i = 0; i < 2; i++) {
            delete [] IEN[i];
            delete [] LM[i];
        }
        delete [] IEN;
        delete [] LM;
        
        // Timing
        t = clock() - t;
        cout << "Previous loop took " << ((float)t/CLOCKS_PER_SEC) << " seconds" << endl;
    }
    
    return 0;
}
