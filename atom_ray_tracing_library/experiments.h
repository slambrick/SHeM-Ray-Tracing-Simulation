/*
 * Copyright (c) 2020, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 */

#ifndef EXPERIMENTS_H_
#define EXPERIMENTS_H_

void generating_rays_cad_pinhole(SourceParam source, int nrays, int *killed,
        int32_t * const cntr_detected, int maxScatters, Surface3D sample, Surface3D plate,
        AnalytSphere the_sphere, double const backWall[], MTRand * const myrng, int32_t * const numScattersRay);

void generating_rays_simple_pinhole(SourceParam source, int n_rays, int * const killed,
        int32_t * const cntr_detected, int maxScatters, Surface3D sample, NBackWall plate,
        AnalytSphere the_sphere, MTRand * const myrng, int32_t * const numScattersRay);

void given_rays_simple_pinhole(Rays3D * const all_rays, int * killed,
        int * const cntr_detected, Surface3D sample, NBackWall plate,
        AnalytSphere the_sphere, int maxScatters, int32_t * const detected,
        int32_t * const which_detector, MTRand * const myrng);

void given_rays_cad_pinhole(Rays3D * const all_rays, int * const killed, int32_t * const cntr_detected,
        Surface3D sample, Surface3D plate, AnalytSphere the_sphere, double const backWall[],
        int maxScatters, int32_t * const detected, MTRand * const myrng);

#endif /* EXPERIMENTS_H_ */
