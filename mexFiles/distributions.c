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
 * This file contains functions that generate new directions for rays upon
 * scattering events. It relies upon the helper functions written in
 * "small_functions.h".
 *
 * Random number generation is performed using the GNU/SL. Functions and code
 * snippets are present that use the C standard library insead.
 */
#include "distributions.h"

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <string.h>

#include "small_functions3D.h"


/*
 * Resolve the distribution function name.
 * This is like an "index" of the scattering distributions,
 * mapping from strings to actual function pointers.
 * If no such name can be found, NULL is returned.
 */
distribution_func distribution_by_name(const char * name) {
    if(strcmp(name, "broad_specular"))
        return broad_specular_scatter;
    if(strcmp(name, "cosine"))
        return cosine_scatter;
    if(strcmp(name, "cosine_specular"))
        return cosine_specular_scatter;
    if(strcmp(name, "uniform"))
        return uniform_scatter;
    return NULL;
}


static double theta_generate(double sigma, gsl_rng *myrng);

static double phi_generate(gsl_rng *myrng);


/*
 * Generate a random direction according to the Gaussian broadened specular.
 *
 * INPUTS:
 *  normal - the surface normal at the ray surface intersection
 *  init_dir - the initial direction of the ray
 *  new_dir  - array to put the new direction in
 *  params   - parameter standard deviation of the scattering distribution
 *  my_rng   - random number genertor object
 */
void broad_specular_scatter(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {

    double theta, phi;
    double theta_normal;
    double sigma = params[1];
    double t0[3];
    double t1[3];
    double t2[3];
    double tester;
    double c_theta;

    /* The 'specular' direction is stored in t0 */
    reflect3D(normal, init_dir, t0);
    /* t1 and t2 are the tangential directions */
    perpendicular_plane(t0, t1, t2);

    do {
        int k;
        double dotted;

        /* Generate a random theta and phi */
        theta = theta_generate(sigma, my_rng);
        phi = phi_generate(my_rng);

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
        tester = gsl_rng_uniform(my_rng);
    } while (c_theta < tester);
    /* We have successfully generated a new direction */
}

/*
 * Generate the polar angle theta according to the Gaussian broaded specular,
 * this angle is the polar angle relative to the specular direction.
 */
static double theta_generate(double sigma, gsl_rng *my_rng) {
    double theta;
    double s_theta = 0;
    double tester = 0;

    /* Sample the Gaussian distribution then reject with a probability that
     * is proportional to sin(theta) */
    do {
        /* Generate a random Gaussian number, note that theta must by non-negative */
        theta = fabs(gsl_ran_gaussian(my_rng, sigma));

        /* If theta exceeds pi then try again */
        if (theta > M_PI)
            continue;

        /* Calculate the sine of the angle */
        s_theta = 0.5*sin(theta);

        /* Generate a tester variable */
        tester = gsl_rng_uniform(my_rng);
    } while (s_theta < tester);

    return(theta);
}

/* Generate a random phi angle for use in the broad specular distribution */
static double phi_generate(gsl_rng *my_rng) {
    double phi;

    phi = 2*M_PI*gsl_rng_uniform(my_rng);
    return(phi);
}

/*
 * Generates a random normalized direction according to a cosine distribution
 * about the provided normal and stores the result in the provided array.
 *
 * INPUTS:
 *  normal  - double array, the normal to the surface at the point of scattering
 *  initial_dir - can be NULL
 *  new_dir - double array, an array to store the new direction of the ray
 *  params  - no parameters expected, can be NULL
 *  my_rng   - gsl_rng pointer, pointer to a GSL random number generator that has
 *            been created and set up with setupGSL()
 */
void cosine_scatter(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {
    double s_theta, c_theta, phi;
    double t1[3];
    double t2[3];

    perpendicular_plane(normal, t1, t2);

    phi = 2*M_PI*gsl_rng_uniform(my_rng);
    s_theta = sqrt(gsl_rng_uniform(my_rng));
    c_theta = sqrt(1 - s_theta*s_theta);

    /* Create the new random direction from the two random angles */
    for (int k = 0; k < 3; k++) {
        new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta + normal[k]*c_theta;
    }
}

/*
 * Generated a random normalized direction according to a cosine distribution
 * about the specular direction and stores the results in the provided array.
 */
void cosine_specular_scatter(const double normal[3], const double initial_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {
    double s_theta, c_theta, phi;
    double theta_normal;
    double t0[3];
    double t1[3];
    double t2[3];


    /* The 'specular' direction is stored in t0 */
    reflect3D(normal, initial_dir, t0);
    perpendicular_plane(t0, t1, t2);

    /* Keep generating direction until one is in the allowed range (not going
     * into the surface */
    do {
        /* Generate random numbers for phi and cos(theta) using the GSL */
        phi = 2*M_PI*gsl_rng_uniform(my_rng);
        s_theta = sqrt(gsl_rng_uniform(my_rng));
        c_theta = sqrt(1 - s_theta*s_theta);

        /* Create the new random direction from the two random angles */
        for (int k = 0; k < 3; k++) {
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
 *  initial_dir - can be NULL
 *  new_dir - double array, an array to store the new direction of the ray
 *  params  - no parameters expected, can be NULL
 *  my_rng  - gsl_rng pointer, pointer to a GSL random number generator that has
 *            been created and set up with setupGSL()
 */
void uniform_scatter(const double normal[3], const double initial_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {
    double s_theta, c_theta, phi;
    double t1[3];
    double t2[3];

    perpendicular_plane(normal, t1, t2);

    /* Generate random numbers for phi and cos(theta) */
    /* Using the GSL */
    phi = 2*M_PI*gsl_rng_uniform(my_rng);
    c_theta = fabs(0.99*gsl_rng_uniform(my_rng) - 1);
    s_theta = sqrt(1 - c_theta*c_theta);

    /* Create the new random direction from the two random angles */
    for (int k = 0; k < 3; k++) {
        new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta + normal[k]*c_theta;
    }
}
