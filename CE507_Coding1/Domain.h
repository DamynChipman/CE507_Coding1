//   Created by Damyn Chipman on 10/24/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   Domain.h
//   PROJECT:   CE507_Coding1

#ifndef Domain_h
#define Domain_h

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

#endif /* Domain_h */
