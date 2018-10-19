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
#include "Linear.h"
#include "FE1D.h"

using namespace std;
using namespace Eigen;

// Global variables
float DEL_X; // To be set later with specific discrete domain
int N[4] = {10,100,1000,10000}; // Number of nodes

// ----- Function Definitions -----
/**
 * @function kElement
 * @brief k Element Matrix
 * @param a : int : Function Index
 * @param b : int : Node Index
 * @returns mat(a,b) : float : Entry a index a,b
 */
float kElement(int a, int b) {
    Eigen::Matrix2f mat;
    
    float value = 1/DEL_X;
    mat(0,0) = value;
    mat(0,1) = -value;
    mat(1,0) = -value;
    mat(1,1) = value;
    
    return mat(a,b);
}

/**
 * @function fElement1
 * @brief f Element Vector for f = c
 * @param a : int : Function Index
 * @param e : int : Element Index
 * @returns vec(a) : float : Entry at index a
 */
float fElement1(int a, int e) {
    Eigen::Vector2f vec;
    float c = 1;
    float f1 = c;
    float f2 = c;
    
    float value = DEL_X/6;
    vec(0) = value*(2*f1 + f2);
    vec(1) = value*(f1 + 2*f2);
    
    return vec(a);
}

/**
 * @function fElement2
 * @brief f Element Vector for f = x
 * @param a : int : Function Index
 * @param e : int : Element Index
 * @returns vec(a) : float : Entry at index a
 */
float fElement2(int a, int e) {
    Eigen::Vector2f vec;
    
    float f1 = (e)*DEL_X;
    float f2 = (e+1)*DEL_X;
    float value = DEL_X/6;
    vec(0) = value*(2*f1 + f2);
    vec(1) = value*(f1 + 2*f2);
    
    return vec(a);
}

/**
 * @function fElement3
 * @brief f Element Vector for f = x^2
 * @param a : int : Function Index
 * @param e : int : Element Index
 * @returns vec(a) : float : Entry at index a
 */
float fElement3(int a, int e) {
    Eigen::Vector2f vec;
    float f1 = (e*DEL_X)*(e*DEL_X);
    float f2 = ((e+1)*DEL_X)*((e+1)*DEL_X);
    
    float value = DEL_X/6;
    vec(0) = value*(2*f1 + f2);
    vec(1) = value*(f1 + 2*f2);
    
    return vec(a);
}

/**
 * @function uActual
 * @brief Actual solution to problem
 * @param x : float : Point to evalute
 * @param p : int : Parameter for f
 * @returns value : float : Value of solution
 */
float uActual(float x, int p) {
    float c = 1;
    if (p == 1) { return 0.5*c*(1 - x*x); }
    else if (p == 2) { return (1.0/6.0)*(1 - x*x); }
    else if (p == 3) { return (1.0/12.0)*(1 - x*x*x*x); }
    else { return 0;}
}

int main(int argc, const char * argv[]) {
    
    // Output file names and N nodes arrays
    string xNames[4] = {"xN1.csv","xN2.csv","xN3.csv","xN4.csv"};
    string u1Names[4] = {"u1N1.csv","u1N2.csv","u1N3.csv","u1N4.csv"};
    string u2Names[4] = {"u2N1.csv","u2N2.csv","u2N3.csv","u2N4.csv"};
    string u3Names[4] = {"u3N1.csv","u3N2.csv","u3N3.csv","u3N4.csv"};
    string uAct1Names[4] = {"uAct1N1.csv","uAct1N2.csv","uAct1N3.csv","uAct1N4.csv"};
    string uAct2Names[4] = {"uAct2N1.csv","uAct2N2.csv","uAct2N3.csv","uAct2N4.csv"};
    string uAct3Names[4] = {"uAct3N1.csv","uAct3N2.csv","uAct3N3.csv","uAct3N4.csv"};
    string IDNames[4] = {"IDN1.csv","IDN2.csv","IDN3.csv","IDN4.csv"};
    string IENNames[4] = {"IENN1.csv","IENN2.csv","IENN3.csv","IENN4.csv"};
    string LMNames[4] = {"LMN1.csv","LMN2.csv","LMN3.csv","LMN4.csv"};
    
    // --- BEGIN LOOP ---
    cout << " --- BEGINNING PROGRAM LOOP --- " << endl;
    clock_t t;   // Timing
    t = clock(); // Timing
    for (int n = 0; n < 4; n++) {
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

        // Perform FE1D
        cout << "   calculating coefs1... " << endl;
        Eigen::VectorXf coefs1 = FE1D(ID, IEN, LM, kElement, fElement1, N[n]); // For f = c
        Eigen::VectorXf coefs2 = FE1D(ID, IEN, LM, kElement, fElement2, N[n]); // For f = x
        Eigen::VectorXf coefs3 = FE1D(ID, IEN, LM, kElement, fElement3, N[n]); // For f = x^2
        
        // Generate solution
        cout << "   generating solutions... " << endl;
        float* u_h1 = (float*)malloc(NPlot * sizeof(float));
        float* u_h2 = (float*)malloc(NPlot * sizeof(float));
        float* u_h3 = (float*)malloc(NPlot * sizeof(float));
        float* u_act1 = (float*)malloc(NPlot * sizeof(float));
        float* u_act2 = (float*)malloc(NPlot * sizeof(float));
        float* u_act3 = (float*)malloc(NPlot * sizeof(float));
        for (int i = 0; i < NPlot; i++) {
            u_h1[i] = 0;
            u_h2[i] = 0;
            u_h3[i] = 0;
            u_act1[i] = uActual(plotDomain.x_h_[i], 1);
            u_act2[i] = uActual(plotDomain.x_h_[i], 2);
            u_act3[i] = uActual(plotDomain.x_h_[i], 3);
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
        ofstream uAct1(uAct1Names[n]);
        ofstream uAct2(uAct2Names[n]);
        ofstream uAct3(uAct3Names[n]);
        ofstream IDOutput(IDNames[n]);
        ofstream IENOutput(IENNames[n]);
        ofstream LMOutput(LMNames[n]);
        char delim = ',';
        
        for (int i = 0; i < NPlot; i++) {
            x << plotDomain.x_h_[i] << delim;
            u1 << u_h1[i] << delim;
            u2 << u_h2[i] << delim;
            u3 << u_h3[i] << delim;
            uAct1 << u_act1[i] << delim;
            uAct2 << u_act2[i] << delim;
            uAct3 << u_act3[i] << delim;
        }
        
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
        free(u_h1);
        free(u_h2);
        free(u_h3);
        free(u_act1);
        free(u_act2);
        free(u_act3);
        
        // Timing
        t = clock() - t;
        cout << "Previous loop took " << ((float)t/CLOCKS_PER_SEC) << " seconds" << endl;
    }
    
    return 0;
}
