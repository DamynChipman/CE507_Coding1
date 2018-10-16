//   Created by Damyn Chipman on 10/6/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   FE1D.h
//   PROJECT:   CE507_Coding1

#ifndef FE1D_h
#define FE1D_h

#include <Eigen/Dense>

//using namespace Eigen;

Eigen::Vector3f FE1D(int* ID, int** IEN, int** LM, float (*k_ab_e)(int,int), float (*f_a_e)(int), int Ne) {
    
    // Initialize matrices
    Eigen::Matrix3f K;
    Eigen::Vector3f F;
    
    // Necessary variables
    float k_ab = 0;
    float f_a = 0;
    int P;
    int Q;
    
    // ----- Begin algorthim -----
    // TODO: Add BC step
    for (int e = 0; e < Ne; e++) {
        for (int a = 0; a < 1; a++) {
            for (int b = 0; b < 1; b++) {
                k_ab = k_ab_e(a,b); // Calc k_ab_e
            }
            f_a = f_a_e(a); // Calc f_a_e
        }

        for (int a = 0; a < 1; a++) {
            P = LM[a][e]; // Integer mapping
            if (P != 0) {
                for (int b = 0; b < 1; b++) {
                    Q = LM[b][e]; // Integer mapping
                    if (Q != 0) {
                        K(P,Q) = K(P,Q) + k_ab;
                        // Update K matrix
                    }
                }
                F(P) = F(P) + f_a;
                // Update F vector
            }
        }
    }

    // Solve system using Eigen
    Eigen::Vector3f d = K.colPivHouseholderQr().solve(F);
    return d;
}

#endif /* FE1D_h */
