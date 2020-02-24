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

void gaussian_specular_scatter(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng);

void broad_specular_scatter(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng);

void cosine_scatter(const double normal[3], const double init_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng);

void cosine_specular_scatter(const double normal[3], const double initial_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng);

void uniform_scatter(const double normal[3], const double initial_dir[3],
        double new_dir[3], const double * params, gsl_rng *my_rng);

#endif
