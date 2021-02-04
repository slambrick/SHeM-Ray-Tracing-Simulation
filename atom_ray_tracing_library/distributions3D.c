/*
 * Copyright (c) 2018-20, Sam Lambrick.
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
#include "distributions3D.h"
#include <stdlib.h>
#include <string.h>
#include "mtwister.h"
#include "common_helpers.h"
#include <math.h>
#include "ray_tracing_core3D.h"

static double theta_generate(double sigma, MTRand * const myrng);

/*
 * Resolve the distribution function name.
 * This is like an "index" of the scattering distributions,
 * mapping from strings to actual function pointers.
 * If no such name can be found, NULL is returned.
 */
distribution_func distribution_by_name(const char * name) {
    if(strcmp(name, "broad_specular") == 0)
        return(diffuse_and_specular);
    if(strcmp(name, "cosine") == 0)
        return(cosine_scatter);
    if(strcmp(name, "cosine_specular") == 0)
        return(cosine_specular_scatter);
    if(strcmp(name, "uniform") == 0)
        return(uniform_scatter);
    if(strcmp(name, "diffraction") == 0)
        return(diffuse_and_diffraction);
    if(strcmp(name, "dw_specular") == 0)
        return(debye_waller_specular);
    if(strcmp(name, "dw_diffraction") == 0)
        return(debye_waller_diffraction);
    if(strcmp(name, "pure_specular") == 0)
        return(pure_specular);
    return NULL;
} 

void pure_specular(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {
    //printf("\nIt has reflected\n");
    reflect3D(normal, init_dir, new_dir);
}

/*
 * Generate rays with broadened specular distribution and a diffuse background.
 *
 * const params:
 *  first the level (0 - 1) of the diffuse background, then sigma of
 * broad_specular.
 */
void diffuse_and_specular(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {

    double diffuse_lvl = params[0];
    double tester;
    genRand(myrng, &tester);
    if(tester < diffuse_lvl)
        cosine_scatter(normal, init_dir, new_dir, params+1, myrng);
    else
        broad_specular_scatter(normal, init_dir, new_dir, params+1, myrng);
}

/*
 * Generate rays according to a 2D diffraction pattern but with cosine-distributed
 * diffuse background.
 *
 * const params:
 *  first the level (0 - 1) of the diffuse background, then as for
 * diffraction_pattern.
 */
void diffuse_and_diffraction(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {

    double diffuse_lvl = params[0];
    double tester;
    genRand(myrng, &tester);
    if(tester < diffuse_lvl)
        cosine_scatter(normal, init_dir, new_dir, params+1, myrng);
    else
        diffraction_pattern(normal, init_dir, new_dir, params+1, myrng);
}


/*
 * Generate rays with some original distribution, and the apply a Debye-Waller
 * factor to the resulting rays. The rays that are rejected are added to the
 * diffuse background.
 *
 * const params:
 *  incident energy in meV
 *  lattice atomic mass (in multiples of proton mass)
 *  lattice temperature in kelvin
 *  lattice Debye temperature in kelvin
 *  std dev of final/initial energy ratio
 * + followed by all the params for the original distribution
 */
void debye_waller_filter_diffuse(distribution_func original_distr,
        const double normal[3], const double init_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {

    // this prefactor appears in the DW exponent if the following
    // are to be in the units stated in the comment above
    const double prefactor = 278.5085;
    double inc_energy = params[0];
    double latt_mass = params[1];
    double temp = params[2];
    double debye_temp = params[3];
    double energy_sigma = params[4];
    double energy_ratio, dwf, tester;

    double exponent = prefactor * inc_energy * temp / latt_mass / (debye_temp * debye_temp);
    gaussian_random_tail(1, energy_sigma, -1, myrng, &energy_ratio);

    // generate a new direction with the original distribution
    original_distr(normal, init_dir, new_dir, params+5, myrng);

    // with probability proportional to debye-waller factor turn it into diffuse scattering
    double tmp;
    dot(init_dir, new_dir, &tmp);
    dwf = exp(- exponent * (1.0 + energy_ratio - 2 * sqrt(energy_ratio) * tmp));
    genRand(myrng, &tester);

    if(tester > dwf)
        cosine_scatter(normal, init_dir, new_dir, NULL, myrng);
}

void debye_waller_specular(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {
    debye_waller_filter_diffuse(broad_specular_scatter, normal, init_dir,
        new_dir, params, myrng);
}


void debye_waller_diffraction(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {
    debye_waller_filter_diffuse(diffraction_pattern, normal, init_dir,
        new_dir, params, myrng);
}

/*
 * Generate rays according to a 2D diffraction pattern given by two
 * reciprocal lattice basis vectors. The general principle is that the
 * incoming (ki) and final (kf) projections of the wave-vectors onto
 * the surface of the sample satisfy kf = ki + g, where g is a
 * linear combination of the two basis vectors.
 *
 * const params:
 *  maximum orders in p and q
 *  a coefficient to pre-multiply the basis vectors
 *  4 floats for 2 x 2D basis vectors
 *  the sigma to broaden the peaks by, and the sigma of the overall gaussian envelope
 */
void diffraction_pattern(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {

    double e1[3], e2[3];    // unit vectors spanning the surface
    double ni[3], nf[3];    // initial and final directions relative to surface
    double delta[2];        // perturbation to smudge the peaks
    double plane_component2;// squared length of the projection parallel to surface
    double tester, gaussian_value; // tester for distributions
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
    dot(init_dir, e1, &ni[0]);
    dot(init_dir, e2, &ni[1]);
    dot(init_dir, normal, &ni[2]);

    do {
        do {
            // generate a random reciprocal vector
            // p is between [-maxp, +maxp], and same for q and maxq
            gen_random_int(2*maxp+1, myrng, &p);
            p = p - maxp;
            gen_random_int(2*maxq+1, myrng, &q);
            q = q - maxq;
            
            // reject to give a Gaussian probability of peaks
            gaussian_value = exp(-(p*p + q*q) / 2 / (envelope_sig*envelope_sig));
            genRand(myrng, &tester);
        } while(tester > gaussian_value);

        // generate gaussian-distributed random perturbation to smudge the peaks
        gaussian_random(0, peak_sig, delta, myrng);

        // add it to the in-plane components of incident direction
        nf[0] = ni[0] + ratio * (p*b1[0] + q*b2[0]) + delta[0];
        nf[1] = ni[1] + ratio * (p*b1[1] + q*b2[1]) + delta[1];
        plane_component2 = nf[0]*nf[0] + nf[1]*nf[1];

    } while(plane_component2 > 1);

    // find the normal component to normalise nf
    nf[2] = sqrt(1 - plane_component2);

    // transform back to lab frame
    for(int i = 0; i < 3; i++)
        new_dir[i] = nf[0] * e1[i] + nf[1] * e2[i] + nf[2] * normal[i];
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
 *  myrng    - 
 */
void broad_specular_scatter(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {

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
        theta = theta_generate(sigma, myrng); // TODO: pointerize
        double uni_rand;
        genRand(myrng, &uni_rand);
        phi = 2*M_PI*uni_rand;

        /* Generate the new direction */
        for (int k = 0; k < 3; k++) {
            new_dir[k] = t1[k]*cos(phi)*sin(theta) + t2[k]*sin(phi)*sin(theta) +
                t0[k]*cos(theta);
        }
        normalise(new_dir);

        /* Calculate the polar angle (normal, new_dir) to the surface normal for
         * the new direction */
        dot(normal, new_dir, &cos_normal);

        /* If the value of theta_normal is greater than pi/2 reject */
        if (cos_normal < 0)
            continue;

        /* We reject the direction with a probability proportional to
         * cos(theta_normal)
         */
        genRand(myrng, &tester);
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
static double theta_generate(double sigma, MTRand * const myrng) {
    double theta;
    double s_theta = 0;
    double tester = 0;
    int cnt;
    double Z[2];
    double rand1;

    /* Sample the Gaussian distribution then reject with a probability that
     * is proportional to sin(theta) */
    cnt = 0;
    do {
        /* Generate a random Gaussian number, note that theta must by
         * non-negative. The Box-muller generates 2 random numbers, store both 
         * and use the previously generated one if we are on an even iteration 
         * (2nd, 4th, etc.)
         */
        if (cnt % 2) {
            gaussian_random(0, sigma, Z, myrng);
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
        genRand(myrng, &tester);
        
        cnt++;
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
 *  myrng   - 
 */
void cosine_scatter(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {
    double s_theta, c_theta, phi;
    double t1[3];
    double t2[3];

    perpendicular_plane(normal, t1, t2);

    double uni_rand;
    genRand(myrng, &uni_rand);
    phi = 2*M_PI*uni_rand;
    genRand(myrng, &uni_rand);
    s_theta = sqrt(uni_rand);
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
        double new_dir[3], const double * const params, MTRand * const myrng) {
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
        /* Generate random numbers for phi and cos(theta) */
    	double uni_rand;
        genRand(myrng, &uni_rand);
        phi = 2*M_PI*uni_rand;
        genRand(myrng, &uni_rand);
        s_theta = sqrt(uni_rand);
        c_theta = sqrt(1 - s_theta*s_theta);

        /* Create the new random direction from the two random angles */
        for (int k = 0; k < 3; k++) {
            new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta +
                t0[k]*c_theta;
        }

        normalise(new_dir);

        /* Calculate the polar angle to the surface normal for the new direction */
        dot(normal, new_dir, &dot_normal);
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
 *  params  - Single number, the maximum allowed polar angle.
 *  myrng   - 
 */
void uniform_scatter(const double normal[3], const double initial_dir[3],
        double new_dir[3], const double * const params, MTRand * const myrng) {
    double s_theta, c_theta, phi;
    double t1[3];
    double t2[3];

    perpendicular_plane(normal, t1, t2);
    double c_theta_min = cos(params[0]);

    /* Generate random numbers for phi and cos(theta) */
    double uni_rand;
    genRand(myrng, &uni_rand);
    phi = 2*M_PI*uni_rand;
    genRand(myrng, &uni_rand);
    c_theta = fabs(c_theta_min*uni_rand - 1);
    s_theta = sqrt(1 - c_theta*c_theta);

    /* Create the new random direction from the two random angles */
    for (int k = 0; k < 3; k++) {
        new_dir[k] = t1[k]*cos(phi)*s_theta + t2[k]*sin(phi)*s_theta + normal[k]*c_theta;
    }
}
