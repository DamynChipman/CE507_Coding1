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
    Eigen::MatrixXf K(Ne,Ne);
    Eigen::VectorXf F(Ne);
    for (int i = 0; i < Ne; i++) {
        for (int j = 0; j < Ne; j++) {
            K(i,j) = 0;
        }
        F(i) = 0;
    }
    
    // Necessary variables
    float k_ab = 0;
    float f_a = 0;
    int P = 0;
    int Q = 0;
    
    std::cout << "--- BEGINNING FE1D LOOP ---" << std::endl;
    // ----- Begin algorthim -----
    // TODO: Add BC step
    for (int e = 0; e < Ne; e++) {
        std::cout << "  e: " << e;
        for (int a = 0; a < 2; a++) {
            std::cout << "  a: " << a;
            for (int b = 0; b < 2; b++) {
                std::cout << "  b: " << b;
                k_ab = k_ab_e(a,b); // Calc k_ab_e
                std::cout << "  k_ab: " << k_ab;
            }
            f_a = f_a_e(a); // Calc f_a_e
            std::cout << "  f_a: " << f_a;
        }
        std::cout << std::endl << "NEXT SECTION" << std::endl;
        for (int a = 0; a < 2; a++) {
            std::cout << "  a: " << a;
            P = LM[a][e]; // Integer mapping
            std::cout << "  P: " << P;
            if (P != -1) {
                for (int b = 0; b < 2; b++) {
                    std::cout << "  b: " << b;
                    Q = LM[b][e]; // Integer mapping
                    std::cout << "  Q: " << Q;
                    if (Q != -1) {
                        K(P,Q) = K(P,Q) + k_ab_e(a,b);
                        std::cout << "  K[P,Q]: " << K(P,Q);
                        // Update K matrix
                    }
                }
                F(P) = F(P) + f_a_e(a);
                std::cout << "  F[P]: " << F(P);
                // Update F vector
            }
        }
        std::cout << std::endl << "END OF LOOP FOR e: " << e << std::endl;
    }
    std::cout << "--- END OF FE1D LOOP ---" << std::endl;
    
    std::cout << K << std::endl;
    std::cout << F << std::endl;

    // Solve system using Eigen
    Eigen::VectorXf d = K.colPivHouseholderQr().solve(F);
    return d;
}

#endif /* FE1D_h */
