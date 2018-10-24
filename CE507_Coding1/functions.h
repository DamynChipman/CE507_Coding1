//   Created by Damyn Chipman on 10/24/18
//   github: @camperD
//   Copyright © 2018 Damyn Chipman. All rights reserved.
//      FILE:   functions.h
//   PROJECT:   CE507_Coding1

#ifndef functions_h
#define functions_h

// Global Variables
float DEL_X; // To be set later with specific discrete domain, needs global scope

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
    //cout << "  In fElement1: e = " << e << endl;
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
    //cout << "  In fElement2: e = " << e << endl;
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
    //cout << "  In fElement3: e = " << e << endl;
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
    else if (p == 2) { return (1.0/6.0)*(1 - x*x*x); }
    else if (p == 3) { return (1.0/12.0)*(1 - x*x*x*x); }
    else { return 0;}
}

#endif /* functions_h */
