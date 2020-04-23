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
 * This fule contains functions that generate new directions for rays upon
 * scattering events. It relies upon the helper functions written in
 * "small_functions.h".
 * 
 * Random number generation is performed using the GNU/SL. Functions and code
 * snippets are present that use the C standard library insead.
 */
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_math.h>
#include "small_functions3D.h"
#include "common_helpers.h"
#include "scattering_processes3D.h"
#include "ray_tracing_structs3D.h"
#include <stdlib.h>
#include <math.h>

static double theta_generate(double sigma);

static double phi_generate();

/* 
 * Generates a new direction for a ray given the initial direction, unit normal 
 * to the surface, and the type of scattering off the sample.
 *
 * INPUTS:
 *  normal     - the surface normal at the ray-surface intersection, 3 element
 *               double array
 *  init_dir   - initial direction of the ray, 3 element double array
 *  new_dir    - 3 element array to put the new direction in
 *  scattering - double describing the scattering off the sample
 *  sigma      - standard deviation of the broad specular distribution, may be
 *               unused.
 */
void new_direction3D(Ray3D *the_ray, double normal[3], double scattering, 
                     double parameters) {
    
    double tester;
    int j;
    double new_dir[3];
    
    
    /* If diffuseLvl is 2 then use a uniform distribution, if it is [0,1] 
     * then use either a cosine or specular distribution, or a combination 
     * of both. If it is 3 then it is cosine scattering about the specular
     */
    if (fabs(scattering - 2) < 0.000001) {
        /* Uniform scattering selected */
        uniformScatter3D(normal, new_dir);
    } else if (fabs(scattering - 3) < 0.000001) {
        /* Specuarly centred cosine scattering selected */
        cosineSpecularScatter3D(normal, the_ray->direction, new_dir);
    } else if (fabs(scattering - 4) < 0.000001) {
        /* Broad Gaussian specular scattering selected */
        broadSpecular3D(normal, the_ray->direction, new_dir, parameters);
    } else {
        /* A combination of diffse cosine and pure specular scattering selected */
        tester = (double)rand() / (double)RAND_MAX;
        if (tester < scattering) {
            cosineScatter3D(normal, new_dir);
        } else {
            reflect3D(normal, the_ray->direction, new_dir);
        }
    }
    
    /* Update the direction of the ray */
    for (j = 0; j < 3; j++)
        the_ray->direction[j] = new_dir[j];
}

/*
 * Generate a random direction according to the Gaussian broadened specular.
 *
 * INPUTS:
 *  normal - the surface normal at the ray surface intersection
 *  init_dir - the initial direction of the ray
 *  new_dir  - array to put the new direction in
 *  sigma    - parameter standard deviation of the scattering distribution
 */
void broadSpecular3D(double normal[3], double init_dir[3], double new_dir[3],
        double sigma) {
    
    double theta, phi;
    double theta_normal;
    double a, b, c;
    double t0[3];
    double t1[3];
    double t2[3];
    double epsilon;
    double tester;
    double c_theta;
    
    epsilon = 1e-3;
    
    /* The 'specular' direction is stored in t0 */
    reflect3D(normal, init_dir, t0);
    
    /* Generate the first tangent vector t1 by requiring t0 dot t1 = 0 
     * and normalise it */
    if (fabs(t0[2]) > epsilon) {
        t1[0] = 1;
        t1[1] = 1; 
        t1[2] = - (t0[0] + t0[1])/t0[2];
    } else if (fabs(t0[1]) > epsilon) {
        t1[0] = 1;
        t1[1] = - (t0[0] + t0[2])/t0[1];
        t1[2] = 1;
    } else {
        t1[0] = - (t0[1] + t0[2])/t0[0];
        t1[1] = 1; 
        t1[2] = 1;
    }
    
    normalise(t1);
    
    /* Generate the second tangent vector t2 using the cross product */
    a = t0[1]*t1[2] - t0[2]*t1[1];
    b = - t0[0]*t1[2] + t0[2]*t1[0];
    c = t0[0]*t1[1] - t0[1]*t1[0];
    t2[0] = a;
    t2[1] = b;
    t2[2] = c;
    
    do {
        int k;
        double dotted;
        
        /* Generate a random theta and phi */
        theta = theta_generate(sigma);
        phi = phi_generate();
        
        /* Generate the new direction */
        for (k = 0; k < 3; k++) {
            new_dir[k] = t1[k]*cos(phi)*sin(theta) + t2[k]*sin(phi)*sin(theta) +
                t0[k]*cos(theta);
        }
        
        normalise(new_dir);
        
        /* Calculate the polar angledouble normal[3], double init_dir[3] to the surface normal for the new direction */
        dotted = dot(normal, new_dir);
        theta_normal = acos(dotted);
        
        /* If the value of theta_normal is greater than pi/2 reject */
        if ((theta_normal >= M_PI/2) || (dotted < 0)) {
            c_theta = 0;
            tester = 1;
            continue;
        }
        
        /* We reject the direction with a probability proportional to 
         * cos(theta_normal) 
         */
        c_theta = cos(theta_normal);
        tester = (double)rand() / (double)RAND_MAX;
    } while (c_theta < tester);
    /* We have successfully generated a new direction */
}

/* 
 * Generate the polar angle theta according to the Gaussian broaded specular,
 * this angle is the polar angle relative to the specular direction.
 */
static double theta_generate(double sigma) {
    double theta;
    double s_theta;
    double tester;
    int cnt;
    double Z[2];
    
    tester = 0;
    s_theta = 0;
    gaussian_random(0, sigma, Z);
    
    /* Sample the Gaussian distribution then reject with a probability that
     * is proportional to sin(theta) */
    cnt = 0;
    do {
        double rand1;
                
        /* Generate a random Gaussian number, note that theta must by
         * non-negative. The Box-muller generates 2 random numbers, store both 
         * and use the previously generated one if we are on an even iteration 
         * (2nd, 4th, etc.)
         */
        if (cnt % 2) {
            gaussian_random(0, sigma, Z);
            rand1 = Z[0];
        } else {
            rand1 = Z[1];
        }
        theta = fabs(rand1);
        
        /* If theta exceeds pi then try again */
        if (theta > M_PI)
            continue;
        
        /* Calculate the sine of the angle */
        s_theta = 0.5*sin(theta);
        
        /* Generate a tester variable */
        tester = (double)rand() / (double)RAND_MAX;
        
        cnt += 1;
    } while (s_theta < tester);
    
    return(theta);
}

/* Generate a random phi angle for use in the broad specular distribution */
static double phi_generate() {
    double phi;
    double rand1;
    
    rand1 = (double)rand() / (double)RAND_MAX;
    phi = 2*M_PI*rand1;
    return(phi);
}

/*
 * Generates a random normalized direction according to a cosine distribution
 * about the provided normal and stores the result in the provided array.
 * 
 * INPUTS:
 *  normal  - double array, the normal to the surface at the point of scattering
 *  new_dir - double array, an array to store the new direction of the ray
 */
void cosineScatter3D(double normal[3], double new_dir[3]) {
    double s_theta, c_theta, phi;
    double a, b, c;
    double t1[3];
    double t2[3];
    int k;
    double epsilon;
    double rand1, rand2;
    
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
    
    rand1 = (double)rand() / (double)RAND_MAX;
    rand2 = (double)rand() / (double)RAND_MAX;
    
    phi = 2*M_PI*rand1;
    s_theta = sqrt(rand2);
    c_theta = sqrt(1 - s_theta*s_theta);
    
    /* Create the new random direction from the two random angles */
    for (k = 0; k < 3; k++) {
        new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta + normal[k]*c_theta;
    }
}

/* 
 * Generated a random normalized direction according to a cosine distribution
 * about the specular direction and stores the results in the provided array.
 */
void cosineSpecularScatter3D(double normal[3], double initial_dir[3], double new_dir[3]) {
    double s_theta, c_theta, phi;
    double theta_normal;
    double a, b, c;
    double t0[3];
    double t1[3];
    double t2[3];
    int k;
    double epsilon;
    double rand1, rand2;
    
    epsilon = 1e-3;
    
    /* The 'specular' direction is stored in t0 */
    reflect3D(normal, initial_dir, t0);
    
    /* Generate the first tangent vector t1 by requiring t0 dot t1 = 0 
     * and normalise it */
    if (fabs(t0[2]) > epsilon) {
        t1[0] = 1;
        t1[1] = 1; 
        t1[2] = - (t0[0] + t0[1])/t0[2];
    } else if (fabs(t0[1]) > epsilon) {
        t1[0] = 1;
        t1[1] = - (t0[0] + t0[2])/t0[1];
        t1[2] = 1;
    } else {
        t1[0] = - (t0[1] + t0[2])/t0[0];
        t1[1] = 1; 
        t1[2] = 1;
    }
    
    normalise(t1);
    
    /* Generate the second tangent vector t2 using the cross product */
    a = t0[1]*t1[2] - t0[2]*t1[1];
    b = - t0[0]*t1[2] + t0[2]*t1[0];
    c = t0[0]*t1[1] - t0[1]*t1[0];
    t2[0] = a;
    t2[1] = b;
    t2[2] = c;
    
    /* Keep generating direction until one is in the allowed range (not going 
     * into the surface */
    do {
        /* Generate random numbers for phi and cos(theta) */
        rand1 = (double)rand() / (double)RAND_MAX;
        rand2 = (double)rand() / (double)RAND_MAX;
        
        phi = 2*M_PI*rand1;
        s_theta = sqrt(rand2);
        c_theta = sqrt(1 - s_theta*s_theta);
        
        /* Create the new random direction from the two random angles */
        for (k = 0; k < 3; k++) {
            new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta +
                t0[k]*c_theta;
        }
        
        normalise(new_dir);
        
        /* Calculate the polar angle to the surface normal for the new direction */
        theta_normal = acos(dot(normal, new_dir));
    } while (theta_normal > M_PI/2);
}

/*
 * Generates a random normalized direction according to a uniform distribution
 * about the provided normal and stores the result in the provided array.
 * 
 * INPUTS:
 *  normal  - double array, the normal to the surface at the point of scattering
 *  new_dir - double array, an array to store the new direction of the ray
 */
void uniformScatter3D(double normal[3], double new_dir[3]) {
    double s_theta, c_theta, phi;
    double a, b, c;
    double t1[3];
    double t2[3];
    int k;
    double epsilon;
    double rand1, rand2;
    
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
    rand1 = (double)rand() / (double)RAND_MAX;
    rand2 = (double)rand() / (double)RAND_MAX;
    phi = 2*M_PI*rand1;
    c_theta = fabs(0.999*rand2 - 1);
    s_theta = sqrt(1 - c_theta*c_theta);
    
    /* Create the new random direction from the two random angles */
    for (k = 0; k < 3; k++) {
        new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta + normal[k]*c_theta;
    }
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
void reflect3D(double normal[3], double init_dir[3], double new_dir[3]) {
    propogate(init_dir, normal, -2*dot(normal, init_dir), new_dir);
    normalise(new_dir);
}
