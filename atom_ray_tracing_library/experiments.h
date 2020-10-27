/*
 * Copyright (c) 2020, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 */

#ifndef EXPERIMENTS_H_
#define EXPERIMENTS_H_

void generating_rays_cad_pinhole(SourceParam source, int nrays, int *killed,
        int32_t * const cntr_detected, int maxScatters, Surface3D const * const sample, Surface3D const * const plate,
        AnalytSphere the_sphere, BackWall * backWall, MTRand * const myrng, int32_t * const numScattersRay);

void generating_rays_simple_pinhole(SourceParam source, int n_rays, int * const killed,
        int32_t * const cntr_detected, int maxScatters, Surface3D const * const sample, NBackWall plate,
        AnalytSphere the_sphere, MTRand * const myrng, int32_t * const numScattersRay);

#endif /* EXPERIMENTS_H_ */
