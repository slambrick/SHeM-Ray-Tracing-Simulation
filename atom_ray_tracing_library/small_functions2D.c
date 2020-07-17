/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 * 
 * Contains a series of small helper functions for performing the scattering of
 * ryas off of a surface in 2D. 
 */
#include "mex.h"
#include "ray_tracing_structs2D.h"
#include "common_helpers.h"
#include "small_functions2D.h"

/* 
 * Calculates the dot product of two two element vectors passed as arrays.
 * 
 * INPUTS:
 *  a - double array, first vector
 *  b - double array, second vector
 *  
 * OUTPUTS:
 *  result - double, the value of the dot product
 */
double dot2(double a[2], double b[2]) {
    int i;
    double result = 0;
    for (i = 0; i < 2; i++) {
        result += a[i]*b[i];
    }
    return result;
}

/* 
 * Returns the square of the norm of a two vector.
 * 
 * INPUTS: 
 *  vect - double array, three element array to calculate the norm^2 of.
 * 
 * OUTPUTS:
 *  result - the value of the norm^2
 */
double norm2(double vect[2]) {
    double result = 0;
    int i;
    for (i = 0; i < 2; i++) {
        result += vect[i]*vect[i];
    }
    return result;
}

/* 
 * Normalizes a 2 element vector. Stores the normalized vector in the same 
 * array as the original.
 * 
 * INPUTS:
 *  vect - double array, three element array to normalize
 */
void normalise2(double vect[2]) {
    double magnitude;
    int i;
    
    magnitude = norm2(vect);
    for (i = 0; i < 2; i++) {
        vect[i] = vect[i]/sqrt(magnitude);
    }
}

/* 
 * Prints out the vertices and normals of the sample surface to the terminal.
 * 
 * INPUTS:
 *  sample - a surface struct containing the vertices and normals of a surface.
 */
void print_surface2D(Surface2D sample) {
    int i;
    
    for (i = 0; i < sample.n_elements; i++) {
        double v1[2];
        double v2[2];
        double nn[2];
        
        /* Get the vertices and normal for the ith element of surface */
        get_element2D(sample, i, v1, v2, nn);
        mexPrintf("Elemnt %i:\n", i + 1);
        mexPrintf("Vertex 1 = ");
        print1D_double(v1, 2);
        mexPrintf("Vertex 2 = ");
        print1D_double(v2, 2);
        mexPrintf("Normal = ");
        print1D_double(nn, 2);
        mexPrintf("\n");
    }
}

/* 
 * Gets the nth element of a Surface struct and puts the information in the 
 * provided arrays.
 * 
 * INPUTS:
 *  Sample - a Surface struct with all the information on the surface in it
 *  n      - integer, the index of the surface element
 *  v1     - two element array for putting the first vertex in
 *  v2     - two element array for putting the second vertex in
 *  nn     - two element array for putting the unit normal in
 */
void get_element2D(Surface2D Sample, int n, double v1[2], double v2[2], 
                   double nn[2]) {
    int m, j;
    
    /* Extract the information about this element of surface */
    m = 2*n;
    for (j = 0; j < 2; j++) {
        v1[j] = Sample.V[m + j];
        v2[j] = Sample.V[m + 2 + j];
        nn[j] = Sample.N[m + j];
    }
}

/* 
 * Puts the rays defined by the arrays of positions and directions into an array
 * of Ray structs.
 * 
 * INPUTS:
 *  rays    - pointer to an array of Ray2D structs
 *  ray_pos - 2 by n array of initial ray positions
 *  ray_dir - 2 by n array of initial ray directions
 *  nray    - the number of rays
 */
void compse_rays2D(Ray2D *rays, double ray_pos[], double ray_dir[], int nrays) {
    int i;
    
    for (i = 0; i < nrays; i++) {
        int n;
        
        n = i*2;
        rays[i].position[0] = ray_pos[n];
        rays[i].position[1] = ray_pos[n+1];
        rays[i].direction[0] = ray_dir[0];
        rays[i].direction[1] = ray_dir[1];
        rays[i].nscatters = 0;
        
        /* Initialise ray not to be on any surface element */
        rays[i].on_element = -1;
    }
    return;
}

/* 
 * Print the position and direction of a ray that is kept inside a ray struct.
 * 
 * INPUT:
 *  ray - a Ray2D struct
 */
void print_ray2D(Ray2D ray) {
    mexPrintf("Position = [%f, %f]\n", ray.position[0], ray.position[1]);
    mexPrintf("Direction = [%f, %f]\n\n", ray.direction[0], ray.direction[1]);
}

/* 
 * Prints the positions and directions of all the rays a 'Rays' object.
 * 
 * INPUT:
 *  all_rays - a Rays object containing a number of rays.
 */
void print_all_rays2D(Rays2D all_rays) {
    int i;
    
    for (i = 0; i < all_rays.nrays; i++) {
        mexPrintf("Ray %i:\n", i + 1);
        print_ray2D(all_rays.rays[i]);
    }
}

/*
 * Puts the directions of all the rays into a single matrix for outputting.
 * 
 * INPUTS:
 *  all_rays  - a Rays object containing a number of rays
 *  final_dir - a two dimensional array of doubles for containing the final
 *              directions of the rays, must be of size 2*nrays where nrays is
 *              the number of rays stored in all_rays
 */
void get_directions2D(Rays2D all_rays, double final_dir[]) {
    int i;
    
    for (i = 0; i < all_rays.nrays; i++) {
        int n;
        
        n = i*2;
        final_dir[n] = all_rays.rays[i].direction[0];
        final_dir[n + 1] = all_rays.rays[i].direction[1];
    }
}

/*
 * Puts the positions of all the rays into a single matrix for outputting.
 * 
 * INPUTS:
 *  all_rays  - a Rays object containing a number of rays
 *  final_dir - a two dimensional array of doubles for containing the final
 *              positions of the rays, must be of size 2*nrays where nrays is
 *              the number of rays stored in all_rays
 */
void get_positions2D(Rays2D all_rays, double final_pos[]) {
    int i;
    
    for (i = 0; i < all_rays.nrays; i++) {
        int n;
        
        n = i*2;
        final_pos[n] = all_rays.rays[i].position[0];
        final_pos[n + 1] = all_rays.rays[i].position[1];
    }
}
