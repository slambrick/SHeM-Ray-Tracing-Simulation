/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 * Function declarations for the small functions used in the SHeM ray tracing
 * simulation.
 */

#ifndef _distributions
#define _distributions

#include <gsl/gsl_rng.h>

/*
 * This is the TYPE of a distribution function. They take in:
 * - a surface normal
 * - an original direction
 * - a new direction, which will be overwritten
 * - a pointer to a double array of parameters
 * - a GSL random number generator
 */
typedef void (*distribution_func)(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng*);

distribution_func distribution_by_name(const char * name);


/*
 * Generate rays with broadened specular distribution and a diffuse background.
 *
 * PARAMS:
 *  first the level (0 - 1) of the diffuse background, then sigma of
 * broad_specular.
 */
void diffuse_and_specular(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng);

/*
 * Generate rays according to a 2D diffraction pattern but with cosine-distributed
 * diffuse background.
 *
 * PARAMS:
 *  first the level (0 - 1) of the diffuse background, then as for
 * diffraction_pattern.
 */
void diffuse_and_diffraction(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng);

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
        double new_dir[3], const double * params, gsl_rng *my_rng);

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
        double new_dir[3], const double * params, gsl_rng *my_rng);

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
        double new_dir[3], const double * params, gsl_rng *my_rng);

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
        double new_dir[3], const double * params, gsl_rng *my_rng);


/*
 * Generated a random normalized direction according to a cosine distribution
 * about the specular direction.
 */
void cosine_specular_scatter(const double normal[3], const double initial_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng);


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
        double new_dir[3], const double * params, gsl_rng *my_rng);

#endif
