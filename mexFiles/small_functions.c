/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Function for using in the ray tracing simulations of the SHeM. Should only
 * be included in files for use with mex and MATLAB.
 * 
 * This file contains a number of small functions that are used in the SHeM ray 
 * tracing simulation, or would be useful in adaptions of that simulation, or 
 * that come in iseful when debugging.
 * 
 * Random number generation is performed using the GNU/SL. Functions and code
 * snippets are present that use the C standard library insead.
 *
 * NOTE: the random number generator seeder must be called before cosieScatter.
 */
#include <gsl/gsl_rng.h>
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Prints out a 3 element 1D double array passed as an argument. */
void print1D(double vect[3]) {
    mexPrintf("{%f, ", vect[0]);
    mexPrintf("%f, ", vect[1]);
    mexPrintf("%f}\n", vect[2]);
}

/* Prints out a 3 element 1D int array passed as an argument. */
void print1Dint(int vect[3]) {
    mexPrintf("{%i, ", vect[0]/3 + 1);
    mexPrintf("%i, ", vect[1]/3 + 1);
    mexPrintf("%i}\n", vect[2]/3 + 1);
}

/* Prints out a 3 by 3 double array passed as an argument. */
void print3x3(double matrix[3][3]) {
    int i;
    mexPrintf("{\n");
    for (i = 0; i < 3; i++) {
        mexPrintf("%f, ", matrix[i][0]);
        mexPrintf("%f, ", matrix[i][1]);
        mexPrintf("%f\n", matrix[i][2]);
    }
    mexPrintf("}\n");
}

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

/* 
 * Reflects a vector through the normal provided storing the result in the 
 * array provided.
 * 
 * INPUTS:
 *  normal   - double array, the normal to the surface at the point of 
 *             reflection
 *  init_dir - double array, the initial direction of the ray
 *  new_dir  - double array, an array to store the final direction of the ray
 */
void reflect(double normal[3], double init_dir[3], double new_dir[3]) {
    propogate(init_dir, normal, -2*dot(normal, init_dir), new_dir);
    normalise(new_dir);
}

/* 
 * Seeds the C standard library random number generator 
 * NOTE: use the GSL random number generator instead.
 */
void seederForScattering(void) {
    srand(time(NULL));
}

/* 
 * Sets up the GNU/SL random number generator and returns a pointer to the
 * generator. This must be called before cosineScatter and the generator this 
 * returns should be passed to cosieScatter.
 * 
 * OUTPUTS:
 *  r - a pointer to a GSL random number generator
 */
gsl_rng* setupGSL(void) {
    gsl_rng *r;
    unsigned long int t;
    
    /*t = (unsigned long int)time(NULL);*/
    
    struct timeval tv;
    gettimeofday(&tv, 0);
    t =  tv.tv_sec + tv.tv_usec;
        
    /* We use the standard random number generator algorithm the Mersenne Twister */
    r = gsl_rng_alloc(gsl_rng_mt19937);
    
    gsl_rng_set(r, t);
    
    return r;
}

/*
 * Generates a random normalized direction according to a cosine distribution
 * about the provided normal and stores the result in the provided array.
 * 
 * INPUTS:
 *  normal  - double array, the normal to the surface at the point of scattering
 *  new_dir - double array, an array to store the new direction of the ray
 *  myrng   - gsl_rng pointer, pointer to a GSL random number generator that has
 *            been created and set up with setupGSL()
 */
void cosineScatter(double normal[3], double new_dir[3], gsl_rng *myrng) {
    double s_theta, c_theta, phi;
    double a, b, c;
    double t1[3];
    double t2[3];
    int k;
    double epsilon;
    
    epsilon = 1e-3;
    
    /* Generate a 'tangent' vector perpendicular to the normal and normalise it */
    if (fabs(normal[2]) > epsilon) {
        t1[0] = 1;
        t1[1] = 1; 
        t1[2] = - (normal[0] + normal[1])/normal[2];
    } else if (fabs(normal[1]) > epsilon) {
        t1[0] = 1;
        t1[1] = - (normal[0] + normal[2])/normal[1];
        t1[2] = 1;
    } else {
        t1[0] = - (normal[1] + normal[2])/normal[0];
        t1[1] = 1; 
        t1[2] = 1;
    }
    
    normalise(t1);
    
    /* Generate the second tangent vector using the cross product */
    a = normal[1]*t1[2] - normal[2]*t1[1];
    b = - normal[0]*t1[2] + normal[2]*t1[0];
    c = normal[0]*t1[1] - normal[1]*t1[0];
    t2[0] = a;
    t2[1] = b;
    t2[2] = c;
    
    /* Generate random numbers for phi and cos(theta) */
    /* Using the standard library */
    if (0) {
        phi = 2*M_PI*((double)rand() / (double)RAND_MAX);
        s_theta = sqrt((double)rand() / (double)RAND_MAX);
        c_theta = sqrt(1 - s_theta*s_theta);
    }
    
    /* Using the GSL */
    if (1) {
        phi = 2*M_PI*gsl_rng_uniform(myrng);
        s_theta = sqrt(gsl_rng_uniform(myrng));
        c_theta = sqrt(1 - s_theta*s_theta);
    }
    
    /* Create the new random direction from the two random angles */
    for (k = 0; k < 3; k++) {
        new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta + normal[k]*c_theta;
    }
}

/*
 * Generates a random normalized direction according to a uniform distribution
 * about the provided normal and stores the result in the provided array.
 * 
 * INPUTS:
 *  normal  - double array, the normal to the surface at the point of scattering
 *  new_dir - double array, an array to store the new direction of the ray
 *  myrng   - gsl_rng pointer, pointer to a GSL random number generator that has
 *            been created and set up with setupGSL()
 */
void uniformScatter(double normal[3], double new_dir[3], gsl_rng *myrng) {
    double s_theta, c_theta, phi;
    double a, b, c;
    double t1[3];
    double t2[3];
    int k;
    double epsilon;
    
    epsilon = 1e-3;
    
    /* Generate a 'tangent' vector perpendicular to the normal and normalise it */
    if (fabs(normal[2]) > epsilon) {
        t1[0] = 1;
        t1[1] = 1; 
        t1[2] = - (normal[0] + normal[1])/normal[2];
    } else if (fabs(normal[1]) > epsilon) {
        t1[0] = 1;
        t1[1] = - (normal[0] + normal[2])/normal[1];
        t1[2] = 1;
    } else {
        t1[0] = - (normal[1] + normal[2])/normal[0];
        t1[1] = 1; 
        t1[2] = 1;
    }
    
    normalise(t1);
    
    /* Generate the second tangent vector using the cross product */
    a = normal[1]*t1[2] - normal[2]*t1[1];
    b = - normal[0]*t1[2] + normal[2]*t1[0];
    c = normal[0]*t1[1] - normal[1]*t1[0];
    t2[0] = a;
    t2[1] = b;
    t2[2] = c;
    
    /* Generate random numbers for phi and cos(theta) */
    /* Using the standard library */
    if (0) {
        phi = 2*M_PI*((double)rand() / (double)RAND_MAX);
        c_theta = fabs(0.99*((double)rand() / (double)RAND_MAX) - 1);
        s_theta = sqrt(1 - c_theta*c_theta);
    }
    
    /* Generate random numbers for phi and cos(theta) */    
    /* Using the GSL */
    if (1) {
        phi = 2*M_PI*gsl_rng_uniform(myrng);
        c_theta = fabs(0.99*gsl_rng_uniform(myrng) - 1);
        s_theta = sqrt(1 - c_theta*c_theta);
    }
    
    /* Create the new random direction from the two random angles */
    for (k = 0; k < 3; k++) {
        new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta + normal[k]*c_theta;
    }
}
