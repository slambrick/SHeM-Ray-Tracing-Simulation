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
    if(strcmp(name, "gaussian_specular") == 0)
        return gaussian_specular_scatter;
    if(strcmp(name, "broad_specular") == 0)
        return diffuse_and_specular;
    if(strcmp(name, "cosine") == 0)
        return cosine_scatter;
    if(strcmp(name, "cosine_specular") == 0)
        return cosine_specular_scatter;
    if(strcmp(name, "uniform") == 0)
        return uniform_scatter;
    if(strcmp(name, "diffraction") == 0)
        return diffuse_and_diffraction;
    return NULL;
}


static double theta_generate(double sigma, gsl_rng *myrng);

/*
 * Generate rays with broadened specular distribution and a diffuse background.
 *
 * PARAMS:
 *  first the level (0 - 1) of the diffuse background, then sigma of
 * broad_specular.
 */
void diffuse_and_specular(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {

    double diffuse_lvl = params[0];
    double tester = gsl_rng_uniform(my_rng);
    if(tester < diffuse_lvl)
        cosine_scatter(normal, init_dir, new_dir, params+1, my_rng);
    else
        broad_specular_scatter(normal, init_dir, new_dir, params+1, my_rng);
}

/*
 * Generate rays according to a 2D diffraction pattern but with cosine-distributed
 * diffuse background.
 *
 * PARAMS:
 *  first the level (0 - 1) of the diffuse background, then as for
 * diffraction_pattern.
 */
void diffuse_and_diffraction(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {

    double diffuse_lvl = params[0];
    double tester = gsl_rng_uniform(my_rng);
    if(tester < diffuse_lvl)
        cosine_scatter(normal, init_dir, new_dir, params+1, my_rng);
    else
        diffraction_pattern(normal, init_dir, new_dir, params+1, my_rng);
}


/*
 * Generate rays according to a 2D diffraction pattern given by two
 * reciprocal lattice basis vectors. The general principle is that the
 * incoming (ki) and final (kf) projections of the wave-vectors onto
 * the surface of the sample satisfy kf = ki + g, where g is a
 * linear combination of the two basis vectors.
 *
 * PARAMS:
 *  maximum orders in p and q
 *  a coefficient to pre-multiply the basis vectors
 *  4 floats for 2 x 2D basis vectors
 *  the sigma to broaden the peaks by, and the sigma of the overall gaussian envelope
 */
void diffraction_pattern(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {

    double e1[3], e2[3];    // unit vectors spanning the surface
    double t1[3], t2[3];    // vectors orthogonal to diffraction peak
    double ni[3], nf[3];    // initial and final directions relative to surface
    double peak[3];         // direction of some diffraction peak
    double theta, phi;      // angles with respect to diffraction peak
    double plane_component2;// squared length of the projection parallel to surface
    double tester, exp_val; // tester for distributions
    int p, q;               // diffraction peak indices

    // unpack the arguments
    const int maxp = (int)params[0], maxq = (int)params[1];
    const double ratio = params[2]; // the lambda/a ratio that scales the reciprocal vector

    double b1[2], b2[2];    // the reciprocal lattice vectors
    b1[0] = params[3];
    b1[1] = params[4];
    b2[0] = params[5];
    b2[1] = params[6];

    double peak_sig = params[7];    // width of individual peaks
    double envelope_sig = params[8]; // width of overall envelope

    // switch to surface-specific coordinates: (x, y) in the plane, z orthogonal:
    perpendicular_plane(normal, e1, e2);
    ni[0] = dot(init_dir, e1);
    ni[1] = dot(init_dir, e2);
    ni[2] = dot(init_dir, normal);

    do {
        do {
            // generate a random reciprocal vector
            // p is between [-maxp, +maxp], and same for q and maxq
            p = gsl_rng_uniform_int(my_rng, 2*maxp+1) - maxp;
            q = gsl_rng_uniform_int(my_rng, 2*maxq+1) - maxq;

            // reject to give a Gaussian probability of peaks
            exp_val = exp(-(p*p + q*q) / 2 / (envelope_sig*envelope_sig));
            tester = gsl_rng_uniform(my_rng);
        } while(tester > exp_val);

        // add it to the in-plane components of incident direction
        peak[0] = ni[0] + ratio * (p*b1[0] + q*b2[0]);
        peak[1] = ni[1] + ratio * (p*b1[1] + q*b2[1]);
        plane_component2 = peak[0]*peak[0] + peak[1]*peak[1];

    } while(plane_component2 > 1);

    // find the normal component to normalise peak
    peak[2] = sqrt(1 - plane_component2);

    // to smudge the peaks, find tangential vectors to this direction
    // and then use the broad_specular algorithm.
    perpendicular_plane(peak, t1, t2);

    do {
        theta = theta_generate(peak_sig, my_rng);
        phi = 2*M_PI*gsl_rng_uniform(my_rng);

        for(int j = 0; j < 3; j++)
            nf[j] = t1[j]*cos(phi)*sin(theta) + t2[j]*sin(phi)*sin(theta) +
                peak[j]*cos(theta);
        normalise(nf);
    // reject the directions into the surface
    // NB in the surface coordinates, z is up
    } while(nf[2] <= 0);

    // transform back to lab frame
    for(int i = 0; i < 3; i++)
        new_dir[i] = nf[0] * e1[i] + nf[1] * e2[i] + nf[2] * normal[i];
}


/*
 * Generate rays with a gaussian-distributed polar angle relative to the specular direction.
 * NB this is not entirely physically realistic, the purpose is mostly to illustrate the need
 * for extra complexity in the broad specular distribution. For actual specular reflections,
 * use broad_specular_scatter instead.
 *
 * INPUTS:
 *  normal - the surface normal at the ray surface intersection
 *  init_dir - the initial direction of the ray
 *  new_dir  - array to put the new direction in
 *  params   - first element must be standard deviation of gaussian distribution
 *  my_rng   - random number generator object
 */
void gaussian_specular_scatter(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {

    double theta, phi, normal_dot;
    double sigma = params[0];
    double specular[3], t1[3], t2[3];
    // double tester;

    // calculate the specular direction
    reflect3D(normal, init_dir, specular);
    // get the tangent vectors
    perpendicular_plane(specular, t1, t2);


    /* NB: the range of coordinates here is theta in [-pi, pi], phi in [0, pi]
     * rather than the usual [0, pi] and [0, 2pi]. The reason is that we are
     * generating theta via a gaussian distribution centred at zero.
     */

    do {
        phi = M_PI*gsl_rng_uniform(my_rng);

        // generate a zero-centred Gaussian distribution with cut-off at pi
        do {
            theta = gsl_ran_gaussian(my_rng, sigma);
        } while(fabs(theta) > M_PI);

        for (int k = 0; k < 3; k++) {
            new_dir[k] = t1[k]*cos(phi)*sin(theta) + t2[k]*sin(phi)*sin(theta) +
                specular[k]*cos(theta);
        }
        normalise(new_dir);

        normal_dot = dot(normal, new_dir);
        // tester = gsl_rng_uniform(my_rng);
    } while(normal_dot < 0);
}


/*
 * Generate a random direction according to the Gaussian broadened specular:
 *
 * The probability distribution is:
 *
 * P = cos(normal) * sin(theta) * exp(-theta^2/(2*sigma^2)),
 * where cos(normal) is the cosine of the new direction with the normal,
 * and theta is the angle between the new direction and the specular.
 *
 * INPUTS:
 *  normal - the surface normal at the ray surface intersection
 *  init_dir - the initial direction of the ray
 *  new_dir  - array to put the new direction in
 *  params   - first element must be standard deviation of gaussian distribution
 *  my_rng   - random number generator object
 */
void broad_specular_scatter(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {

    double theta, phi;
    double cos_normal;
    double sigma = params[0];
    double t0[3], t1[3], t2[3];
    double tester = 1.0;

    /* The 'specular' direction is stored in t0 */
    reflect3D(normal, init_dir, t0);
    /* t1 and t2 are the tangential directions */
    perpendicular_plane(t0, t1, t2);

    do {
        /* Generate a random theta and phi */
        theta = theta_generate(sigma, my_rng);
        phi = 2*M_PI*gsl_rng_uniform(my_rng);

        /* Generate the new direction */
        for (int k = 0; k < 3; k++) {
            new_dir[k] = t1[k]*cos(phi)*sin(theta) + t2[k]*sin(phi)*sin(theta) +
                t0[k]*cos(theta);
        }
        normalise(new_dir);

        /* Calculate the polar angle (normal, new_dir) to the surface normal for the new direction */
        cos_normal = dot(normal, new_dir);

        /* If the value of theta_normal is greater than pi/2 reject */
        if (cos_normal < 0)
            continue;

        /* We reject the direction with a probability proportional to
         * cos(theta_normal)
         */
        tester = gsl_rng_uniform(my_rng);
    } while (cos_normal < tester);
    /* We have successfully generated a new direction */
}

/*
 * Generate the polar angle theta according to the Gaussian broadened specular,
 * this angle is the polar angle relative to the specular direction.
 *
 * The distribution is P = sin(theta) * exp(-theta^2/(2*sigma^2)),
 * where the sine comes from the solid angle integrated over azimuthal directions.
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
 * about the specular direction.
 */
void cosine_specular_scatter(const double normal[3], const double initial_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng) {
    double s_theta, c_theta, phi;
    double dot_normal;
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
        dot_normal = dot(normal, new_dir);
    } while (dot_normal < 0);
}

/*
 * Generates a random normalized direction with uniformly-distributed polar
 * angle about the provided normal.
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
    c_theta = fabs(0.999*gsl_rng_uniform(my_rng) - 1);
    s_theta = sqrt(1 - c_theta*c_theta);

    /* Create the new random direction from the two random angles */
    for (int k = 0; k < 3; k++) {
        new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta + normal[k]*c_theta;
    }
}