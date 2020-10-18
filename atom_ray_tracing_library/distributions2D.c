/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 * 
 * Contains functions for creating a new direction for a scattered ray. To
 * include a new scattering distribution add the direction generating function
 * here and add a new case to the function new_direction in scattering2D (or
 * in the future scattering3D).
 * 
 * Functions should be of the form:
 *  scatter(normal, init_dir, new_dir, distribution parameters..., random number
 *      generator)
 */
#include "mex.h"
#include "small_functions2D.h"
#include "distributions2D.h"
#include <stdlib.h>
#include "mtwister.h"
#include <math.h>

/* 
 * Gives the ray a new direction given the normal to the surface. The direction
 * of the ray is updated here. There are different ways of updating the
 * direction depending on the scattering distribution desired. Parameters are
 * passed in as an array of doubles, passing parameters as an array allows 
 * distribution to have arbitraty number of parametres.
 *     At present the only scattering distribution that uses a parameter is the
 * broad specular distribution which takes a single parameter.
 * 
 * INPUTS:
 *  the_ray    - a Ray2D struct containing the direction and postion of the ray
 *  normal     - a 2 element double array containing the unit normal to the surface
 *               from which to generate the new direction
 *  scattering - int, the type of scattering to model:
 *               0 - specular reflection
 *               1 - diffuse cosine scattering
 *               2 - completely uniform scattering
 *               3 - broad specular distribution
 * parameters  - double array, parameters for the scattering distribution
 */
void new_direction2D(Ray2D *the_ray, double normal[2], int scattering, 
        double *parameters, MTRand *my_rng) {
    double new_dir[2];
    
    /* Choose the correct type of scattering */
    switch (scattering) {
        case 0:
            reflect2D(normal, the_ray->direction, new_dir);
            break;
        case 1:
            scatterCosine2D(normal, new_dir, my_rng);
            break;
        case 2:
            scatterUniform2D(normal, new_dir, my_rng);
            break;
        case 3:
            scatterBroadSpecular2D(normal, the_ray->direction, new_dir,
                parameters[0], my_rng);
            break;
    }
    
    /* Update the ray struct */
    the_ray->direction[0] = new_dir[0];
    the_ray->direction[1] = new_dir[1];
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
void reflect2D(double normal[2], double init_dir[2], double new_dir[2]) {
    int i;
    double a;
    
    a = -2*dot2(normal, init_dir);
    for (i = 0; i < 2; i++) {
        new_dir[i] = init_dir[i] + normal[i]*a;
    }
    
    /* Ensure that the new direction is normalised. */
    normalise2(new_dir);
}

/*
 * Scatters a ray according to a cosine distribution.
 * 
 * INPUTS: 
 *  normal   - double array, the normal to the surface at the point of 
 *             reflection
 *  init_dir - double array, the initial direction of the ray
 *  new_dir  - double array, an array to store the final direction of the ray
 */
void scatterCosine2D(double normal[2], double new_dir[2], MTRand *my_rng) {
    double s_theta;
    double c_theta;
    double t[2];
    double epsilon;
    
    epsilon = 1e-3;
    
    /* Create a tangent vector. */
    if (fabs(normal[1]) > epsilon) {
        t[0] = 1;
        t[1] = normal[0]*t[0]/normal[1];
    } else {
        t[1] = 1;
        t[0] = normal[1]*t[1]/normal[0];
    }
    normalise2(t);
    
    /* Generate random angle. */
    double tmp;
    genRand(my_rng, &tmp);
    s_theta = 2*tmp - 1;
    c_theta = sqrt(1 - s_theta*s_theta);
    
    /* Create the new direction */
    new_dir[0] = t[0]*s_theta + normal[0]*c_theta;
    new_dir[1] = t[1]*s_theta + normal[1]*c_theta;
}

/*
 * Scatters a ray according to a cosine distribution.
 * 
 * INPUTS: 
 *  normal   - double array, the normal to the surface at the point of 
 *             reflection
 *  init_dir - double array, the initial direction of the ray
 *  new_dir  - double array, an array to store the final direction of the ray
 */
void scatterUniform2D(double normal[2], double new_dir[2], MTRand *my_rng) {
    double theta;
    double s_theta;
    double c_theta;
    double t[2];
    double epsilon;
    
    epsilon = 1e-3;
    
    /* Create a tangent vector. */
    if (fabs(normal[1]) > epsilon) {
        t[0] = 1;
        t[1] = normal[0]*t[0]/normal[1];
    } else {
        t[1] = 1;
        t[0] = normal[1]*t[1]/normal[0];
    }
    normalise2(t);
    
    /* Generate random angle. */
    double tmp;
    genRand(my_rng, &tmp);
    theta = (2*tmp - 1)*M_PI/2;
    s_theta = sin(theta);
    c_theta = cos(theta);
    
    /* Create the new direction */
    new_dir[0] = t[0]*s_theta + normal[0]*c_theta;
    new_dir[1] = t[1]*s_theta + normal[1]*c_theta;
}

/*
 * Scatters a ray according to a cosine distribution.
 * 
 * INPUTS: 
 *  normal   - double array, the normal to the surface at the point of 
 *             reflection
 *  init_dir - double array, the initial direction of the ray
 *  new_dir  - double array, an array to store the final direction of the ray
 *  sigma_P  - double, the 'parameter standard deviation' of the distirbution
 */
void scatterBroadSpecular2D(double normal[2], double init_dir[2], 
        double new_dir[2], double sigma_P, MTRand *my_rng) {
    /* TODO */
}
