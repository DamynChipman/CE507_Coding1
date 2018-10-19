//   Created by Damyn Chipman on 10/6/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   FE1D.h
//   PROJECT:   CE507_Coding1

#ifndef FE1D_h
#define FE1D_h

#include <Eigen/Dense>

//using namespace Eigen;

Eigen::VectorXf FE1D(int* ID, int** IEN, int** LM, float (*k_ab_e)(int,int), float (*f_a_e)(int,int), int Ne) {
    
    // Initialize matrices
    Eigen::MatrixXf K(Ne,Ne);
    Eigen::VectorXf F(Ne);
    for (int i = 0; i < Ne; i++) {
        for (int j = 0; j < Ne; j++) {
            K(i,j) = 0;
        }
        F(i) = 0;
    }
    
    // Necessary variables
    int P = 0;
    int Q = 0;
    
    // ----- Begin algorthim -----
    for (int e = 0; e < Ne; e++) {
        for (int a = 0; a < 2; a++) {
            P = LM[a][e]; // Integer mapping
            if (P != -1) {
                for (int b = 0; b < 2; b++) {
                    Q = LM[b][e]; // Integer mapping
                    if (Q != -1) {
                        K(P,Q) = K(P,Q) + k_ab_e(a,b); // Update K matrix
                    }
                }
                F(P) = F(P) + f_a_e(a,e); // Update F vector
            }
        }
    }
    // Apply BC
    P = LM[0][0];
    F(P) = F(P) + f_a_e(0,0);

    // Solve system using Eigen
    Eigen::VectorXf d(Ne);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXf> dec(K);
    d = dec.solve(F);
    return d;
}

#endif /* FE1D_h */
