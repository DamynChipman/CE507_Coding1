//   Created by Damyn Chipman on 10/6/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   FE1D.h
//   PROJECT:   CE507_Coding1

#ifndef FE1D_h
#define FE1D_h

#include "/Users/Damyn/Documents/Coding Files/XCode/Numerical Data Structures/Numerical Data Structures/Matrix.h"

Matrix<float> FE1D(int* ID, int** IEN, int** LM, float (*k_ab_e)(int,int), float (*f_a_e)(int), int Ne) {
    
    // Initialize matrices
    Matrix<float> d(Ne,1);
    Matrix<float> K(Ne,Ne);
    Matrix<float> F(Ne,1);
    
    // Necessary variables
    float k_ab;
    float f_a;
    int P;
    int Q;
    float toAdd;
    
    // ----- Begin algorthim -----
    // Element iterator
    for (int e = 0; e < Ne; e++) {
        //
        for (int a = 0; a < 1; a++) {
            //
            for (int b = 0; b < 1; b++) {
                k_ab = k_ab_e(a,b); // Calc k_ab_e
            }
            f_a = f_a_e(a); // Calc f_a_e
        }
        
        //
        for (int a = 0; a < 1; a++) {
            P = LM[a][e]; // 
            if (P != 0) {
                for (int b = 0; b < 1; b++) {
                    Q = LM[b][e];
                    if (Q != 0) {
                        toAdd = K.at(P,Q);
                        K.set(P, Q, toAdd + k_ab);
                    }
                }
                toAdd = F.at(P,1);
                F.set(P, 1, toAdd + f_a);
            }
        }
    }
    
    // Solve system for d
    
    return d;
};


#endif /* FE1D_h */
