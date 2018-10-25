//   Created by Damyn Chipman on 10/24/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   gaussQuad.h
//   PROJECT:   CE507_Coding1

#ifndef gaussQuad_h
#define gaussQuad_h

#include <cmath>

float gaussQuad(float* u, float* u_h, int Ne) {
    
    // Variables
    //float var = sqrt(3.0/5.0);
    //float evalPoints[3] = {-var, 0, var};
    float weights[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    float error = 0;
    for (int e = 0; e < Ne; e++) {
        for (int i = 0; i < 3; i++) {
            float A = u[e] - u_h[e];
            float B = abs(A);
            float C = pow(B,2);
            float D = DEL_X/2;
            float E = weights[i];
            float F = pow(C*D*E,0.5);
            float toAdd = F;
            //float toAdd = pow(pow(abs(u[e] - u_h[e]),2)*(DEL_X/2)*weights[i],0.5);
            error = error + toAdd;
        }
    }
    
    return error;
}

#endif /* gaussQuad_h */
