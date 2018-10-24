//   Created by Damyn Chipman on 10/6/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   Linear.h
//   PROJECT:   CE507_Coding1

#ifndef Linear_h
#define Linear_h

#include "BasisFunction.h"
#include "Domain.h"

class Linear : public BasisFunction {
    
public:
    
    class Domain {
    public:
        Domain(float left, float right, int n) : lBound_(left), rBound_(right), N_(n) {
            x_h_ = (float*)malloc(n * sizeof(float));
            
            del_x_ = (right - left)/(n);
            for (int i = 0; i < n; i++) {
                x_h_[i] = i*del_x_;
            }
        }
        
        // Accessor Functions
        float getLBound() { return lBound_; }
        float getRBound() { return rBound_; }
        int getN() { return N_; }
        float* getPoints() { return x_h_; }
        float getDelX() { return del_x_; }
        
        ~Domain() {
            //            free(x_h_);
        }
        
    private:
        float lBound_;
        float rBound_;
        int N_;
        float* x_h_;
        float del_x_;
    };
    
    Linear(float xL, float xR, int N) : domain(Domain(xL,xR,N)) {}
    
    float eval(int A, float x) {
        float xAm1 = domain.getPoints()[A-1];
        float xAp1 = domain.getPoints()[A+1];
        float xA = domain.getPoints()[A];
        float res;
        if (x > xAm1 && x <= xA) {
            res = (x - xAm1)/domain.getDelX();
        }
        else if (x >= xA && x < xAp1) {
            res = (xAp1 - x)/domain.getDelX();
        }
        else {
            res = 0;
        }
        return res;
    }
    
    Domain getDomain() { return domain; }
    
    ~Linear() {}
    
private:
    Domain domain;
};

#endif /* Linear_h */
