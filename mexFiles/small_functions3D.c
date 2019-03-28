/*
 * Copyright (c) 2018-19, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Function for using in the ray tracing simulations of the SHeM. Should only
 * be included in files for use with mex and MATLAB.
 * 
 * This file contains a number of small functions that are used in the SHeM ray 
 * tracing simulation, or would be useful in adaptions of that simulation.
 */
#include "mex.h"
#include "common_helpers.h"
#include "small_functions3D.h"
#include <gsl/gsl_math.h>

/* 
 * Adds a 3 element array to another 3 element array multiplied by a scalar 
 * (propagates an array).
 * 
 * EQUATION:
 *  result = init + direc*a
 * 
 * INPUTS:
 *  init   - double array, vector 1 (initial position)
 *  direc  - double array, vector 2 (direction)
 *  a      - double, amount to multiply vector 2 by (move in the direction of
 *           vector direc by this amount.) 
 *  result - double array, array to store the result (final position) in
 */
void propogate(double init[3], double direc[3], double a, double result[3]) {
    int i;
    for (i = 0; i < 3; i++) {
        result[i] = init[i] + direc[i]*a;
    }
}

/* 
 * Calculates the dot product of two three element vectors passed as arrays.
 * 
 * INPUTS:
 *  a - double array, first vector
 *  b - double array, second vector
 *  
 * OUTPUTS:
 *  result - double, the value of the dot product
 */
double dot(double a[3], double b[3]) {
    int i;
    double result = 0;
    for (i = 0; i < 3; i++) {
        result += a[i]*b[i];
    }
    return result;
}

/* 
 * Returns the square of the norm of a three vector.
 * 
 * INPUTS: 
 *  vect - double array, three element array to calculate the norm^2 of.
 * 
 * OUTPUTS:
 *  result - the value of the norm^2
 */
double norm2(double vect[3]) {
    double result = 0;
    int i;
    for (i = 0; i < 3; i++) {
        result += vect[i]*vect[i];
    }
    return result;
}

/* 
 * Normalizes a 3 element vector. Stores the normalized vector in the same 
 * array as the original.
 * 
 * INPUTS:
 *  vect - double array, three element array to normalize
 */
void normalise(double vect[3]) {
    double magnitude;
    int i;
    
    magnitude = norm2(vect);
    for (i = 0; i < 3; i++) {
        vect[i] = vect[i]/sqrt(magnitude);
    }
}

/* 
 * Solves a 3D matrix equation Au=v for u, we write
 *      (a b c)   (u[0])   (j)
 *      (d e f) . (u[1]) = (k)
 *      (g h i)   (u[2])   (l)
 * Uses Cramer's rule and explicit formula for the determinants to solve the 
 * equation. Checks to see if the determinant of the matrix is less than the 
 * supplied tolerance epsilon, returns 0 if it is less and 1 if it is more than
 * the tolerance.
 * 
 * INPUTS:
 *  A       - double array, the 3x3 matrix, as an array, for the equation
 *  u       - double array, a 3 element array for the answer to be put in
 *  v       - double array, the 3-vector on the other side of the equation
 *  epsilon - double, the tolerance for the size of the determinant of the 
 *            matrix A
 * 
 * OUTPUTS:
 *  result - int, 0 or 1 depending on the size of the determinant of A
 * 
 * NOTE: This function is not used anymore, the equation is explicity solved 
 *       within larger functions for the sake of speed improvments. This
 *       function is kept as it would make any new code written much more 
 *       readable. The commented code here uses a series of variables to make
 *       the function more readable, it is kept for refernce.
 */
int solve3x3(double A[3][3], double u[3], double v[3], double epsilon) {
    double M, Dx, Dy, Dz, X1, X2, X3;
    /*double a,b,c,d,e,f,g,h,i,j,k,l;
    
    a = A[0][0];
    b = A[0][1];
    c = A[0][2];
    
    d = A[1][0];
    e = A[1][1];
    f = A[1][2];
    
    g = A[2][0];
    h = A[2][1];
    i = A[2][2];
    
    j = v[0];
    k = v[1];
    l = v[2];*/
    
    /* These are used twice - only compute them once */
    /*X1 = e*i - h*f;
    X2 = c*h - b*i;
    X3 = b*f - c*e;
    M = a*X1 + d*X2 + g*X3;*/
    
    X1 = A[1][1]*A[2][2] - A[2][1]*A[1][2];
    X2 = A[0][2]*A[2][1] - A[0][1]*A[2][2];
    X3 = A[0][1]*A[1][2] - A[0][2]*A[1][1];
    M = A[0][0]*X1 + A[1][0]*X2 + A[2][0]*X3;
    
    /* fabs() is the math.h abs function for floats */
    if (fabs(M) < epsilon) {
        return 0;
    }
    /*Dx = j*X1 + k*X2 + l*X3;
    
    X1 = a*k - j*d;
    X2 = j*g - a*l;
    X3 = d*l - k*g;
    
    Dy = i*X1 + f*X2 + c*X3;
    Dz = - h*X1 - e*X2 - b*X3;*/
    
    Dx = v[0]*X1 + v[1]*X2 + v[2]*X3;
    
    X1 = A[0][0]*v[1] - v[0]*A[1][0];
    X2 = v[0]*A[2][0] - A[0][0]*v[2];
    X3 = A[1][0]*v[2] - v[1]*A[2][0];
    
    Dy = A[2][2]*X1 + A[1][2]*X2 + A[0][2]*X3;
    Dz = - A[2][1]*X1 - A[1][1]*X2 - A[0][1]*X3;
    
    u[0] = Dx/M;
    u[1] = Dy/M;
    u[2] = Dz/M;
    
    return 1;
}
