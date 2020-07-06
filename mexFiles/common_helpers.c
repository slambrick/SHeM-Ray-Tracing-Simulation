/*
 * Copyright (c) 2018-19, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 * Contains small helper functions common to both 2D and 3D ray Tracing.
 */
#include "common_helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mtwister.h"

/*
 * Linearise [row][column] coordinates in an array of coordinates
 * such as a 3-column n-row matrix of vertices of a Surface3D
 */
int lin(int row, int col) {
    return 3*row + col;
}

/*
 * Prints a 2 or 3 element 1D array to the terminal.
 *
 * INPUTS:
 *  vect - a double pointer to either a 2 or 3 element vector
 *  dim  - int, the dimension of the vector, 2 or 3
 */
void print1D_double(double *vect, int dim) {
    if (dim == 2)
        printf("[%f, %f]\n", vect[0], vect[1]);
    else if (dim == 3)
        printf("[%f, %f, %f]\n", vect[0], vect[1], vect[2]);
    else
        printf("dim into print1D must be 2 or 3.");
}

/*
 * Prints a 2 or 3 element 1D array to the terminal.
 *
 * INPUTS:
 *  vect - an int pointer to either a 2 or 3 element vector
 *  dim  - int, the dimension of the vector, 2 or 3
 */
void print1D_int(int *vect, int dim) {
    if (dim == 2)
        printf("[%i, %i]\n", vect[0], vect[1]);
    else if (dim == 3)
        printf("[%i, %i, %i]\n", vect[0], vect[1], vect[2]);
    else
        printf("dim into print1D must be 2 or 3.");
}

/* Prints out a 3 by 3 double array passed as an argument. */
void print3x3(double matrix[3][3]) {
    int i;
    printf("{\n");
    for (i = 0; i < 3; i++) {
        printf("%f, ", matrix[i][0]);
        printf("%f, ", matrix[i][1]);
        printf("%f\n", matrix[i][2]);
    }
    printf("}\n");
}

/* 
 * Generates two gaussian random numbers using the box-muller transform.
 */
void gaussian_random(double mu, double sigma, double Z[2], MTRand *myrng) {
    double U1, U2;
    
    U1 = genRand(myrng);
    U2 = genRand(myrng);
    
    Z[0] = sqrt(-2*log(U1))*cos(2*M_PI*U2);
    Z[1] = sqrt(-2*log(U1))*sin(2*M_PI*U2);
    
    Z[0] = Z[0]*sigma + mu;
    Z[1] = Z[1]*sigma + mu;
    
    return;
}

/* 
 * Samples from the top tail of a Gaussian distribution. Samples from the
 * distribution and then checks to see if the value is below the cutoff.
 */
double gaussian_random_tail(double mu, double sigma, double cutoff, MTRand *myrng) {
    double Z[2];
    double rand1 = 0;
    int cnt = 0;
    
    cutoff = -1.0;
    do {
        if (cnt % 2) {
            gaussian_random(mu, sigma, Z, myrng);
            rand1 = Z[0];
        } else {
            rand1 = Z[1];
        }
        
        cnt++;
    } while (rand1 < cutoff);
    
    return rand1;
}

/*
 * Create a random int in the desired range, 0 to max-1
 */
int gen_random_int(int max, MTRand *myrand) {
    double uniform_rand;
    int random_int;
    
    uniform_rand = genRand(myrand);
    uniform_rand = uniform_rand*(double)max;
    random_int = (int)floor(uniform_rand);
    
    return random_int;
}
